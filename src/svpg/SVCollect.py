import logging
import pysam

from svpg.SVSignature import SignatureDeletion, SignatureInsertion
from svpg.util import analyze_cigar_indel, merge_cigar

def decompose_cigars(alignment, bam, query_name, min_length):
    """Parse BAM record to extract SVs from CIGAR."""
    sv_signatures = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    indels = analyze_cigar_indel(alignment.cigartuples, min_length)
    read_seq = alignment.query_sequence
    for pos_ref, pos_read, length, typ in indels:
        start = ref_start + pos_ref
        if typ == "DEL":
            sv_signatures.append(SignatureDeletion(ref_chr, start, length, "cigar", query_name, read_seq=read_seq, pos_read=pos_read))
        elif typ == "INS":
            sv_signatures.append(SignatureInsertion(ref_chr, start, length, "cigar", query_name, read_seq=read_seq, pos_read=pos_read))

    return sv_signatures

def analyze_split_indel(alignment_current, alignment_next, ultra_ins_flag=False):
    """Parse BAM record to extract SVs from inter-alignment."""
    distance_on_read = alignment_next['q_start'] - alignment_current['q_end']
    if not alignment_current['is_reverse']:
        distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
        if alignment_next['is_reverse']:  # INV:+-
            if alignment_current['ref_end'] > alignment_next['ref_end']:
                distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_start']
            else:
                distance_on_reference = alignment_current['ref_end'] - alignment_next['ref_start']
    else:
        distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
        if not alignment_next['is_reverse']:  # INV:-+
            if alignment_current['ref_end'] > alignment_next['ref_end']:
                distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_start']
            else:
                distance_on_reference = alignment_current['ref_end'] - alignment_next['ref_start']

    split_signature = []
    deviation = distance_on_read - distance_on_reference
    if alignment_current['ref_chr'] == alignment_next['ref_chr']:
        if alignment_current['is_reverse'] == alignment_next['is_reverse']:
            if distance_on_reference >= -50:
                # INS
                if deviation >= 50:
                    if not alignment_current['is_reverse']:
                        if not ultra_ins_flag:
                            start = alignment_current['ref_end']
                        else:
                            start = min(alignment_current['ref_end'], alignment_next['ref_start'])
                        pos_read = alignment_current['q_end']
                    else:
                        if not ultra_ins_flag:
                            start = alignment_current['ref_start']
                        else:
                            start = min(alignment_current['ref_start'], alignment_next['ref_end'])
                        pos_read = alignment_next['infer_read_length'] - alignment_next['q_start']
                    split_signature.append(SignatureInsertion(
                        alignment_current['ref_chr'], start, deviation, "suppl", alignment_current['read_name'], read_seq=alignment_current['atgc_seq'], pos_read=pos_read, alt_seq='<INS>'))
                # DEL
                elif deviation <= -50:
                    if not alignment_current['is_reverse']:
                        start = alignment_current['ref_end']
                        pos_read = alignment_current['q_end']
                    else:
                        start = alignment_next['ref_end']
                        pos_read = alignment_next['infer_read_length'] - alignment_next['q_start']
                    split_signature.append(SignatureDeletion(alignment_current['ref_chr'], start, -deviation, "suppl", alignment_current['read_name'], read_seq=alignment_current['atgc_seq'], pos_read=pos_read))
                else:
                    return split_signature

    return split_signature

def decompose_split(primary, supplementaries, bam):
    """Parse BAM record to extract SVs from split_reads."""
    read_name = primary.query_name
    alignments = [primary] + supplementaries
    alignment_list = []
    sig_list = []
    for alignment in alignments:
        if alignment.is_reverse:
            q_start = alignment.infer_read_length() - alignment.query_alignment_end
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        alignment_dict = {
            'read_name': read_name,
            'q_start': q_start,
            'q_end': q_end,
            'ref_chr': bam.getrname(alignment.reference_id),
            'ref_start': alignment.reference_start,
            'ref_end': alignment.reference_end,
            'is_reverse': alignment.is_reverse,
            'mapping_quality': alignment.mapping_quality,
            'infer_read_length': alignment.infer_read_length(),
            'atgc_seq': primary.query_sequence,
        }
        alignment_list.append(alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    for index in range(len(sorted_alignment_list) - 1):
        sig_list.extend(analyze_split_indel(sorted_alignment_list[index], sorted_alignment_list[index + 1]))
    if len(alignment_list) >= 3 and sorted_alignment_list[0]['ref_chr'] != sorted_alignment_list[1]['ref_chr']:
        sig_list.extend(
            analyze_split_indel(sorted_alignment_list[0], sorted_alignment_list[-1], ultra_ins_flag=True))

    return sig_list

def retrieve_other_alignments(main_alignment, bam):
    """Reconstruct other alignments of the same read for a given alignment from the SA tag"""
    if main_alignment.get_cigar_stats()[0][5] > 0:
        return []
    try:
        sa_tag = main_alignment.get_tag("SA").split(";")
    except KeyError:
        return []
    other_alignments = []
    # For each other alignment encoded in the SA tag
    for element in sa_tag:
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = main_alignment.query_name
        a.query_sequence = main_alignment.query_sequence
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        a.cigarstring = cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = main_alignment.query_qualities
        a.set_tags([("NM", nm, "i")])

        other_alignments.append(a)

    return other_alignments

def read_bam(contig, start, end, options):
    """Parse BAM record to extract SVs."""
    bam = pysam.AlignmentFile(options.bam, threads=options.num_threads)
    sv_signatures, sv_signatures_inter = [], []
    for current_alignment in bam.fetch(contig, start, end):
        try:
            if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq or current_alignment.reference_start < start:
                continue
            sigs = decompose_cigars(current_alignment, bam, current_alignment.query_name, options.min_sv_size)
            if sigs:
                sigs = merge_cigar(sigs, options.read)
                sv_signatures.extend(sigs)
            if not current_alignment.is_supplementary:
                supplementary_alignments = retrieve_other_alignments(current_alignment, bam)
                good_suppl_alns = [aln for aln in supplementary_alignments if
                                   not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]
                sig_list = decompose_split(current_alignment, good_suppl_alns, bam)
                sv_signatures_inter.extend(sig_list)

        except StopIteration:
            break
        except KeyboardInterrupt:
            logging.warning('Execution interrupted by user. Stop detection and continue with next step..')
            break

    sv_signatures = sv_signatures + sv_signatures_inter

    return sv_signatures
