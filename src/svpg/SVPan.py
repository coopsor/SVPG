import re
from collections import defaultdict

from svpg.SVSignature import SignatureDeletion, SignatureInsertion
from svpg.util import analyze_cigar_indel, merge_cigar

CIGAR_PATTERN = re.compile(r'(\d+)([MIDNSHP=X])')

class Gaf:
    def __init__(self):
        self.query_name = ""
        self.query_length = 0
        self.query_start = 0
        self.query_end = 0
        self.strand = ""
        self.path = ""
        self.path_length = 0
        self.path_start = 0
        self.path_end = 0
        self.mapping_quality = 0
        self.is_primary = True
        self.cigar = ""
        self.ds = ""

def parse_gaf_line(tokens, gfa_node):
    """Parse a single GAF line"""
    gafline = Gaf()

    if '@' in tokens[0]:
        bam_tags = tokens[0].split('@')
        gafline.query_name = bam_tags[0]
        gafline.type = bam_tags[1]
        gafline.pos = bam_tags[2]
    else:
        gafline.query_name = tokens[0]

    gafline.query_length = int(tokens[1])
    gafline.query_start = int(tokens[2])
    gafline.query_end = int(tokens[3])
    gafline.path = re.findall(r'[<>][^<>]+', tokens[5])
    gafline.strand = '+' if gafline.path[0][0] == '>' else '-'
    # gafline.strand = tokens[4]
    gafline.contig = gfa_node[gafline.path[0][1:]].contig
    gafline.offset = gfa_node[gafline.path[0][1:]].offset
    gafline.path_length = int(tokens[6])
    gafline.path_start = int(tokens[7])
    gafline.path_end = int(tokens[8])
    gafline.mapping_quality = int(tokens[11])
    gafline.is_primary = True

    for tok in tokens:
        if "tp:A:" in tok:
            if tok[5:7] != "P":
                gafline.is_primary = False
        if "cg:Z:" in tok:
            gafline.cigar = tok[5:]
        if "ds:Z:" in tok:
            gafline.ds = tok[6:]

    return gafline

def decompose_split(g_list):
    """Parse GAF record to extract SVs from split_reads."""
    alignment_list, split_signature = [], []
    for g in g_list:
        strand = g.strand
        ref_start = g.offset + g.path_start
        ref_end = g.offset + g.path_end
        q_start = g.query_start
        q_end = g.query_end
        alignment_dict = {
            'read_name': g.query_name,
            'q_start': q_start,
            'q_end': q_end,
            'ref_id': g.contig,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'is_reverse': strand == '-'
        }
        alignment_list.append(alignment_dict)
    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    ultra_ins_flag = True if len(alignment_list) >= 3 and sorted_alignment_list[0]['ref_id'] != sorted_alignment_list[1]['ref_id'] else False

    for alignment_current, alignment_next in zip(sorted_alignment_list[:-1], sorted_alignment_list[1:]):
        distance_on_read = alignment_next['q_start'] - alignment_current['q_end']
        if not alignment_current['is_reverse']:
            distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
            if alignment_next['is_reverse']:  # INV:+-
                if alignment_current['ref_end'] > alignment_next['ref_end']:
                    distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_start']
                else:
                    distance_on_reference = alignment_current['ref_end'] - alignment_next['ref_start']
        else:
            continue
            # distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
            # if not alignment_next['is_reverse']:  # INV:-+
            #     if alignment_current['ref_end'] > alignment_next['ref_end']:
            #         distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_start']
            #     else:
            #         distance_on_reference = alignment_current['ref_end'] - alignment_next['ref_start']
        deviation = distance_on_read - distance_on_reference

        if alignment_current['ref_id'] == alignment_next['ref_id']:
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                if distance_on_reference >= -50:
                    if deviation >= 50:  # INS
                        if not alignment_current['is_reverse']:
                            if not ultra_ins_flag:
                                start = alignment_current['ref_end']
                            else:
                                start = min(alignment_current['ref_end'], alignment_next['ref_start'])
                        else:
                            if not ultra_ins_flag:
                                start = alignment_current['ref_start']
                            else:
                                start = min(alignment_current['ref_start'], alignment_next['ref_end'])
                        split_signature.append(
                            SignatureInsertion(alignment_current['ref_id'], start, deviation, "suppl",
                                               alignment_current['read_name'], alt_seq='<INS>',
                                               qry_start=alignment_current['q_end'], qry_end=alignment_next['q_start']))
                    elif deviation <= -50:  # DEL
                        if not alignment_current['is_reverse']:
                            start = alignment_current['ref_end']
                        else:
                            start = alignment_next['ref_end']
                        split_signature.append(
                            SignatureDeletion(alignment_current['ref_id'], start, -deviation, "suppl", alignment_current['read_name']))
                    else:
                        continue

    return split_signature

def pan_node_offset(pan_node, node_list, gfa_node):
    """ Find contig and coordinates for pan_node according to linear_node """
    pan_len = 0
    for node in node_list[node_list.index(pan_node):]:
        node_id = node[1:]
        if gfa_node[node[1:]].sr == 0:
            node_contig = gfa_node[node_id].contig
            node_offset = gfa_node[node_id].offset - pan_len
            return (node_contig, node_offset)
        pan_len += gfa_node[node_id].len
    else:
        return None

def check_continuity(node_list):
    node_num = [int(s[2:]) for s in node_list][::-1]
    if all(node_num[i] - node_num[i - 1] == 1 for i in range(1, len(node_num))):
        return False
    else:
        return node_list

def extract_tsd_alt(ds_seq):
    # ds_seq: '[aatttttgtattt]ttaa...'
    pattern = r'(?:\[([^\[\]]+)\])?([^\[\]]*)(?:\[([^\[\]]+)\])?'
    match = re.fullmatch(pattern, ds_seq)
    if match:
        tsdl, alt, tsdr = match.groups()
        return tsdl or "", alt or "", tsdr or ""
    else:
        return "", ds_seq, ""

def get_node_index_for_pos(pos, cum_lengths):
    for i in range(len(cum_lengths)-1):
        if cum_lengths[i] <= pos < cum_lengths[i+1]:
            return i
    return None

def decompose_cigars(g, gfa_node, min_indel_length=40):
    sigs = []
    node_list = g.path  # ['>s1','>s2']
    first_node = node_list[0]
    first_node_len = gfa_node[first_node[1:]].len
    if gfa_node[first_node[1:]].sr != 0:
        node_result = pan_node_offset(first_node, node_list, gfa_node)
        if node_result:
            g.contig, g.offset = node_result
        else:  # node list is composed entirely of pan_nodes
            return []

    parsed_cigar = CIGAR_PATTERN.findall(g.cigar)
    cigar_tuple = [(int(length), operation) for length, operation in parsed_cigar]
    vars = analyze_cigar_indel(cigar_tuple, min_indel_length, is_gaf=True)
    if vars:
        DS_PATTERN = re.compile(rf'[+-]([atcgn\[\]]{{{min_indel_length},}})', re.IGNORECASE)
        parsed_ds = DS_PATTERN.findall(g.ds) if g.ds else []
    else:
        return []

    effective_first_len = first_node_len - g.path_start
    cum_lengths = [0, effective_first_len]
    for node in node_list[1:]:
        node_len = gfa_node[node[1:]].len
        cum_lengths.append(cum_lengths[-1] + node_len)

    ref_chr = g.contig
    global_ref = g.offset + g.path_start
    last_found_node_index = 0
    for var_index, (pos_ref, pos_read, length, typ) in enumerate(vars):
        ds_str = parsed_ds[var_index] if parsed_ds else ''
        ltsd, alt_seq, rtsd = extract_tsd_alt(ds_str)
        ltsd_len, alt_len, rtsd_len = len(ltsd), len(alt_seq), len(rtsd)

        if (ltsd or rtsd) and length < 1000:  # TSD
            indel_left = pos_ref - ltsd_len
            indel_right = pos_ref + alt_len + rtsd_len

            left_node_index = get_node_index_for_pos(indel_left, cum_lengths)
            right_node_index = get_node_index_for_pos(indel_right - 1, cum_lengths)

            if left_node_index != right_node_index:  # Filter cross-node indels
                continue

        start = None
        if g.strand == '+':  # first node is forward
            start = global_ref + pos_ref
        else:
            # if check_continuity(node_list):
            #     return []
            if pos_ref < first_node_len - g.path_start:  # the indel is in first node
                global_ref = g.offset + (first_node_len - g.path_start)
                if typ == 'INS':
                    start = global_ref - pos_ref
                else:
                    start = global_ref - pos_ref - length
            else:  # find global_ref according to node offset
                start_index = max(0, last_found_node_index)
                for i in range(start_index, len(node_list)):
                    node = node_list[i]
                    node_name = node[1:]
                    node_len = gfa_node[node_name].len

                    node_start_ref = cum_lengths[i]
                    node_end_ref = cum_lengths[i + 1]

                    if node_start_ref <= pos_ref < node_end_ref:
                        last_found_node_index = i

                        if gfa_node[node_name].sr == 0:  # the node is a linear node
                            global_ref = gfa_node[node_name].offset
                        else:  # the node is a pan node
                            node_result = pan_node_offset(node, node_list, gfa_node)
                            if not node_result:
                                break
                            else:
                                global_ref = node_result[1]

                        local_ref = pos_ref - node_start_ref
                        if node[0] == '>':  # the node is forward
                            start = global_ref + local_ref
                        else:  # the node is reverse
                            if typ == "INS":
                                start = global_ref + node_len - local_ref
                            else:
                                start = global_ref + node_len - local_ref - length
                        break
        if start is None:
            continue

        if typ == "DEL":
            sigs.append(SignatureDeletion(ref_chr, start, length, "cigar", g.query_name))
        elif typ == "INS":
            sigs.append(SignatureInsertion(ref_chr, start, length, "cigar", g.query_name, alt_seq=ltsd+alt_seq+rtsd))

    return sigs


def read_gaf(gaf_chunk, gfa_node, options):
    """Parse SVsignatures GAF record to extract SVs."""
    inconsist_read = []
    sv_signatures = []
    read_dict = {}
    j, k, e = 0, 0, 0
    min_sv_size = options.min_sv_size
    for i, line in enumerate(gaf_chunk):
        tokens = line.strip().split('\t')
        if tokens[4] == '*':
            continue

        g = parse_gaf_line(tokens, gfa_node)
        if g.mapping_quality < options.min_mapq:
            j += 1
            continue

        if tokens[0] in read_dict:
            read_dict[tokens[0]].append(g)
        else:
            read_dict[tokens[0]] = [g]

        node_list = g.path  # ['>s1','>s2','>s3']
        node_sr = [gfa_node[node[1:]].sr for node in node_list]
        sigs = []

        if sum(node_sr) == 0:
            for node_current, node_next in zip(node_list[:-1], node_list[1:]):
                # map to non-adjacent nodes, ['>s1', '>s3']
                split_node_temp = list(range(min(int(node_current[2:]), int(node_next[2:])) + 1,
                                             max(int(node_current[2:]), int(node_next[2:]))))
                if len(split_node_temp) > 0:
                    if node_current[0] == '>':  # ['>s1', '>s3']
                        start = gfa_node[node_current[1:]].offset + gfa_node[node_current[1:]].len
                        end = gfa_node[node_next[1:]].offset
                    else:  # ['<s3', '>s1']
                        start = gfa_node[node_next[1:]].offset
                        end = gfa_node[node_current[1:]].offset
                    sigs.append(SignatureDeletion(g.contig, start, end - start, "ref_split", g.query_name, node_ls=split_node_temp))

            sigs_cigar = decompose_cigars(g, gfa_node, min_sv_size)
            sigs.extend(sigs_cigar)
        else:
            sigs_liner_pan = decompose_cigars(g, gfa_node, min_indel_length=10)
            sigs_cigar = [sig for sig in sigs_liner_pan if sig.svlen >= min_sv_size]

            linear_index = [i for i, x in enumerate(node_sr) if x == 0]
            linear_node = [node_list[i] for i in linear_index]
            for node_current, node_next in zip(linear_node[:-1], linear_node[1:]):
                if node_current[0] == '>':
                    start = gfa_node[node_current[1:]].offset + gfa_node[node_current[1:]].len
                    # Only retaining the cigar SVs of the linear nodes
                    sigs_cigar = [sig for sig in sigs_cigar if
                                  not (start <= sig.start <= gfa_node[node_next[1:]].offset)]
                else:
                    start = gfa_node[node_next[1:]].offset
                    sigs_cigar = [sig for sig in sigs_cigar if
                                  not (start <= sig.start <= gfa_node[node_current[1:]].offset)]
                split_node_temp = list(range(min(int(node_current[2:]), int(node_next[2:])) + 1,
                                             max(int(node_current[2:]), int(node_next[2:]))))
                split_len = sum([gfa_node['s' + str(node)].len for node in split_node_temp])
                pan_node = node_list[node_list.index(node_current) + 1:node_list.index(node_next)]

                # map to a pan node: insersion
                if pan_node:
                    length = sum([gfa_node[pan[1:]].len for pan in pan_node])
                    for indel in sigs_liner_pan:
                        if start <= indel.start <= start + length:  # cigar SVs in the pan node
                            if indel.type == "DEL":
                                length -= indel.svlen
                            else:
                                length += indel.svlen
                    if length - split_len >= min_sv_size:
                        alt_seq = ''.join([gfa_node[node[1:]].sequence for node in pan_node])
                        sigs.append(SignatureInsertion(g.contig, start, length - split_len, "ref_split", g.query_name, alt_seq=alt_seq, node_ls=pan_node))

                if len(split_node_temp) > 0:  # map to a missing linear nodes: deletion
                    if sum([gfa_node['s' + str(node)].len for node in split_node_temp]) < min_sv_size:
                        continue
                    if node_current[0] == '>':
                        end = gfa_node[node_next[1:]].offset
                    else:
                        end = gfa_node[node_current[1:]].offset
                    if pan_node and length >= min_sv_size:  # the pan node inserted in the middle
                        end -= length
                    if end - start >= min_sv_size:
                        sigs.append(SignatureDeletion(g.contig, start, end - start, "ref_split", g.query_name))

            sigs = sigs+sigs_cigar

        sigs_initial = [sig for sig in sigs if sig.type == g.type]
        sigs = merge_cigar(sigs_initial, options.read)

        bam_pos = g.pos.split(':')
        bam_len = int(bam_pos[2]) - int(bam_pos[1])
        sigs_ = []
        # Find the closest SV record
        if len(sigs) > 1:
            min_distance = float("inf")
            for sig in sigs:
                sig_distance = abs(int(bam_pos[1]) - sig.start) + abs(int(bam_pos[2]) - sig.end)
                if sig_distance < min_distance:
                    min_distance = sig_distance
                    sigs_ = [sig]
        elif len(sigs) == 1:
            sigs_ = sigs

        if sigs_ and min(sigs_[0].svlen, bam_len) / max(sigs_[0].svlen, bam_len) < 0.5:  # todo:remap
            if tokens[0].split('@')[-1] == 'merged_indel':
                pass
            else:
                inconsist_read.append(tokens[0])
                k += 1
                continue

        sv_signatures.extend(sigs)

    for key, value in read_dict.items():
        if len(value) > 1:
            e += 1
            var_split = decompose_split(value)
            sv_signatures.extend(var_split)

    print(f"low_mapq: {j}, inconsistent_read: {k},  split_read: {e}")

    return sv_signatures


def read_gaf_pan(gaf_chunk, gfa_node, options):
    """Parse WGS GAF record to extract SVs."""
    sv_signatures = []
    read_dict = defaultdict(list)

    for line in gaf_chunk:
        tokens = line.strip().split('\t')
        if tokens[4] == '*':
            continue
        g = parse_gaf_line(tokens, gfa_node)
        if g.mapping_quality < options.min_mapq:
            continue

        read_dict[tokens[0]].append(g)

        if g.query_end - g.query_start < g.query_length * 0.7:  # filter cigar in short alignments
            continue

        sigs = decompose_cigars(g, gfa_node, options.min_sv_size)
        if len(sigs) > 1 and len(sigs) > g.query_length * 1e-4 * 2:
            continue
        sigs_merged = merge_cigar(sigs, options.read)
        sv_signatures.extend(sigs_merged)

    for key, value in read_dict.items():
        if len(value) > 1:
            var_split = decompose_split(value)
            sv_signatures.extend(var_split)

    return sv_signatures

