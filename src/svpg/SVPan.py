import re
from collections import defaultdict
import numpy as np

from svpg.SVSignature import SignatureDeletion, SignatureInsertion, SignatureDuplicationTandem, SignatureInversion, SignatureTranslocation
from svpg.util import analyze_cigar_indel, merge_cigar, chr_to_sort_key

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
    try:
        gafline.strand = '+' if gafline.path[0][0] == '>' else '-'
        # gafline.strand = tokens[4]
    except IndexError:
        raise ValueError(f"Please check the GAF file format. SVPG expects standard GAF format, refer to readme for GFA and rGFA format.")
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

def decompose_split(g_list, gfa_node):
    """Parse GAF record to extract SVs from split_reads."""
    alignment_list, split_signature = [], []
    g_list = sorted(g_list, key=lambda g: g.query_start)

    # === Filter short alignments at both ends ===
    zen = len(g_list)
    for j in range(len(g_list) - 1, -1, -1):
        y = g_list[j]
        aln_len = y.query_end - y.query_start
        if aln_len < 2000:
            zen = j
        else:
            break
    if zen < 2:
        return []

    zst = 0
    for j in range(zen):
        y = g_list[j]
        aln_len = y.query_end - y.query_start
        if aln_len < 2000:
            zst = j + 1
        else:
            break
    if zen - zst < 2:
        return []

    g_list = g_list[zst:zen]
    for g in g_list:
        if g.contig == 'pan':
            continue
        strand_start = g.strand
        strand_end = '+' if g.path[-1][0] == '>' else '-'

        ref_start = g.offset + g.path_start if strand_start == '+' else g.offset + gfa_node[g.path[0][1:]].len - g.path_start
        # ref_end = g.offset + g.path_end if strand_end == '+' else g.offset + gfa_node[g.path[-1][1:]].len - g.path_end
        if strand_end == '+':
            ref_end = g.offset + g.path_end
        else:
            if len(g.path) == 1:
                ref_end = g.offset + gfa_node[g.path[-1][1:]].len - g.path_end
            else:
                ref_end = gfa_node[g.path[-1][1:]].offset + gfa_node[g.path[-1][1:]].len - (g.path_end - sum(gfa_node[n[1:]].len for n in g.path[:-1]))
        ref_start, ref_end = min(ref_start, ref_end), max(ref_start, ref_end)

        alignment_dict = {
            'read_name': g.query_name,
            'q_start': g.query_start,
            'q_end': g.query_end,
            'ref_id': g.contig,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'is_reverse_start': strand_start == '-',
            'is_reverse_end': strand_end == '-',
        }
        alignment_list.append(alignment_dict)
    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    ultra_ins_flag = True if len(alignment_list) >= 3 and sorted_alignment_list[0]['ref_id'] != sorted_alignment_list[1]['ref_id'] else False

    for alignment_current, alignment_next in zip(sorted_alignment_list[:-1], sorted_alignment_list[1:]):
        distance_on_read = alignment_next['q_start'] - alignment_current['q_end']
        ref_chr = alignment_current['ref_id']
        read_name = alignment_current['read_name']
        if alignment_current['ref_id'] == alignment_next['ref_id']:
            if alignment_current['is_reverse_end'] == alignment_next['is_reverse_start']:
                if not alignment_current['is_reverse_end'] and not alignment_next['is_reverse_start']:  # ++
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
                else:  # --
                    distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
                if distance_on_reference >= -50:
                    deviation = distance_on_read - distance_on_reference
                    if deviation >= 50:  # INS
                        if not alignment_current['is_reverse_end']:
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
                                               ))
                    elif deviation <= -50:  # DEL
                        if not alignment_current['is_reverse_end']:
                            start = alignment_current['ref_end']
                        else:
                            start = alignment_next['ref_end']
                        split_signature.append(
                            SignatureDeletion(alignment_current['ref_id'], start, -deviation, "suppl", alignment_current['read_name']))
                    else:
                        continue
                else:  # DUP
                    # if distance_on_reference <= -options.min_sv_size:
                    if not alignment_current['is_reverse_end']:
                        start = alignment_next['ref_start']
                        end = alignment_current['ref_end']
                    else:
                        start = alignment_current['ref_start']
                        end = alignment_next['ref_end']

                    sv_sig = (alignment_current['ref_id'], start, end, "suppl", alignment_current['read_name'])
                    split_signature.append(SignatureDuplicationTandem(*sv_sig))
            else:  # INV
                mid_c = (alignment_current['ref_start'] + alignment_current['ref_end']) / 2
                c_len = alignment_current['ref_end'] - alignment_current['ref_start']
                if not alignment_current['is_reverse_end']:  # +-
                    overlap = min(alignment_current['ref_end'], alignment_next['ref_end']) - max(
                        alignment_current['ref_start'], alignment_next['ref_start'])

                    # case1: next entirely left
                    if alignment_next['ref_end'] <= alignment_current['ref_start'] - 5:
                        start, end = alignment_next['ref_end'], alignment_current['ref_end']
                        label = 'left_rev'
                        strand1, strand2 = 'rev', 'fwd'
                    # case2: right-half(inv) overlap
                    elif overlap > 0 and alignment_next['ref_end'] > mid_c and overlap / c_len >= 0.5:
                        start, end = min(alignment_next['ref_end'], alignment_current['ref_end']), max(alignment_next['ref_end'], alignment_current['ref_end'])
                        label = 'left_rev_foldback'
                        strand1, strand2 = 'rev', 'fwd'
                    # case3: left-half overlap
                    elif overlap > 0 and alignment_next['ref_start'] < mid_c and overlap / c_len >= 0.5:
                        start, end = min(alignment_next['ref_end'], alignment_current['ref_end']), max(alignment_next['ref_end'], alignment_current['ref_end'])
                        label = 'left_fwd_foldback'
                        strand1, strand2 = 'fwd', 'rev'
                    # case4: next entirely right
                    elif alignment_next['ref_start'] >= alignment_current['ref_end'] + 5:
                        start, end = alignment_current['ref_end'], alignment_next['ref_end']
                        label = 'left_fwd'
                        strand1, strand2 = 'fwd', 'rev'
                    else:
                        continue
                else:
                    overlap = min(alignment_current['ref_end'], alignment_next['ref_end']) - max(
                        alignment_current['ref_start'], alignment_next['ref_start'])
                    # case1: current entirely left
                    if alignment_next['ref_start'] >= alignment_current['ref_end'] + 5:
                        start, end = alignment_current['ref_start'], alignment_next['ref_start']
                        label = 'right_fwd'
                        strand1, strand2 = 'rev', 'fwd'
                    # case2: foldback-right (right half overlap)
                    elif overlap > 0 and alignment_next['ref_end'] > mid_c and overlap / c_len >= 0.5:
                        start, end = min(alignment_next['ref_start'], alignment_current['ref_start']), max(alignment_next['ref_start'], alignment_current['ref_start'])
                        label = 'right_fwd_foldback'
                        strand1, strand2 = 'rev', 'fwd'
                    # case3: foldback-left (left half overlap)
                    elif overlap > 0 and alignment_next['ref_start'] < mid_c and overlap / c_len >= 0.5:
                        start, end = min(alignment_next['ref_start'], alignment_current['ref_start']), max(alignment_next['ref_start'], alignment_current['ref_start'])
                        label = 'right_rev_foldback'
                        strand1, strand2 = 'fwd', 'rev'
                    # case4: current entirely right
                    elif alignment_next['ref_end'] <= alignment_current['ref_start'] - 5:
                        start, end = alignment_next['ref_start'], alignment_current['ref_start']
                        label = 'right_rev'
                        strand1, strand2 = 'fwd', 'rev'
                    else:
                        continue

                svsize = end - start
                if svsize >= 50:
                    sv_sig = (ref_chr, start, end, "suppl", read_name, label)
                    split_signature.append(SignatureInversion(*sv_sig))

        else:  # TRA
            ref_chr_next = alignment_next['ref_id']
            ref_chr_key = chr_to_sort_key(ref_chr)
            ref_chr_next_key = chr_to_sort_key(ref_chr_next)
            if not ref_chr_key or not ref_chr_next_key:
                continue
            if alignment_current['is_reverse_end'] == alignment_next['is_reverse_start']:
                #  (++, --)
                if not alignment_current['is_reverse_end']:  # ++
                    strand1, strand2 = 'fwd', 'fwd'
                    start = alignment_current['ref_end']
                    end = alignment_next['ref_start']
                else:  # --
                    strand1, strand2 = 'rev', 'rev'
                    start = alignment_current['ref_start']
                    end = alignment_next['ref_end']
            else:
                # (+-, -+)
                if not alignment_current['is_reverse_end']:  # +-
                    strand1, strand2 = 'fwd', 'rev'
                    start = alignment_current['ref_end']
                    end = alignment_next['ref_end']
                else:  # -+
                    strand1, strand2 = 'rev', 'fwd'
                    start = alignment_current['ref_start']
                    end = alignment_next['ref_start']

            # --- ensure consistent chr order ---
            if ref_chr_key > ref_chr_next_key:
                ref_chr, ref_chr_next = ref_chr_next, ref_chr
                start, end = end, start
                if (strand1, strand2) == ('fwd', 'fwd'):
                    strand1, strand2 = 'rev', 'rev'
                elif (strand1, strand2) == ('rev', 'rev'):
                    strand1, strand2 = 'fwd', 'fwd'
                else:
                    strand1, strand2 = strand1, strand2

            sv_sig = (ref_chr, start, strand1, ref_chr_next, end, strand2, "suppl", read_name)
            split_signature.append(SignatureTranslocation(*sv_sig))

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

def decompose_cigars(g, gfa_node, options, min_indel_length=50):
    sigs = []
    node_list = g.path  # ['>s1','>s2']

    first_node = node_list[0]
    first_node_len = gfa_node[first_node[1:]].len
    if gfa_node[first_node[1:]].sr != 0:
        g.contig = 'pan'
        return []

    hap_contigs = set()
    for node in node_list:
        node_contig = gfa_node[node[1:]].contig
        hap_contigs.add(node_contig)
    if options.read == 'hifi':
        if g.strand == '-' and len(hap_contigs) >= 2:## if g.strand == '-':return []
            return []
    else:
        if len(hap_contigs) >= 3:
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

def calculate_euclidean_distance_sigs(sig1, sig2, weights=[1, 5]):
    """Calculate Euclidean distance between two SV signatures."""
    sig1 = np.array(sig1, dtype=float)
    sig2 = np.array(sig2, dtype=float)

    diff = sig2 - sig1
    dist = np.sqrt(np.sum(weights * diff ** 2, axis=1))

    return dist


def read_gaf(gaf_chunk, gfa_node, options):
    """Parse SVsignatures GAF record to extract SVs."""
    sv_signatures = []
    read_dict = {}
    j, k, e = 0, 0, 0
    min_sv_size = options.min_sv_size
    for i, line in enumerate(gaf_chunk):
        tokens = line.strip().split('\t')
        if tokens[4] == '*':
            continue

        g = parse_gaf_line(tokens, gfa_node)
        bam_pos = g.pos.split(':')
        bam_len = int(bam_pos[2]) - int(bam_pos[1])

        if g.mapping_quality < options.min_mapq:
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
                    sigs.append(SignatureDeletion(g.contig, start, end - start, "ref_split", g.query_name, pan_node=split_node_temp))

            sigs_cigar = decompose_cigars(g, gfa_node, options, min_sv_size)
            sigs.extend(sigs_cigar)
        else:
            sigs_liner_pan = decompose_cigars(g, gfa_node, options, min_indel_length=10)
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

                # map to a pan node: insertion
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
                        sigs.append(SignatureInsertion(g.contig, start, length - split_len, "ref_split", g.query_name, alt_seq=alt_seq, pan_node=pan_node))

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

        sigs = [sig for sig in sigs if sig.type == g.type]

        sigs_ = []
        # Find the closest SV record
        if len(sigs) > 1:
            dis = calculate_euclidean_distance_sigs([bam_pos[1], bam_len], [[sig.start, sig.svlen] for sig in sigs])
            min_index = np.argmin(dis)
            sigs_ = [sigs[min_index]]
        elif len(sigs) == 1:
            sigs_ = sigs

        if not sigs_ or (sigs_ and min(sigs_[0].svlen, bam_len) / max(sigs_[0].svlen, bam_len) < 0.7):
            if g.type == "INS":
                sv_signatures.append(SignatureInsertion(bam_pos[0], int(bam_pos[1]), bam_len, "inconsistent", g.query_name, alt_seq="<INS>"))
            else:
                sv_signatures.append(SignatureDeletion(bam_pos[0], int(bam_pos[1]), bam_len, "inconsistent", g.query_name))
            continue

        sv_signatures.extend(sigs_)

    for key, value in read_dict.items():
        if len(value) > 1:
            e += 1
            var_split = decompose_split(value, gfa_node)
            sv_signatures.extend(var_split)

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

        sigs = decompose_cigars(g, gfa_node, options)
        if len(sigs) > 1 and len(sigs) > g.query_length * 1e-4 * 2:
            continue
        sigs_merged = merge_cigar(sigs, max_merge=options.max_merge_threshold)
        sv_signatures.extend(sigs_merged)

    for key, value in read_dict.items():
        if len(value) > 1:
            var_split = decompose_split(value, gfa_node)
            sv_signatures.extend(var_split)

    return sv_signatures