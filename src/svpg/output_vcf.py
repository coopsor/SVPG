import numpy as np
import time
import re
import os.path
from collections import defaultdict, Counter

import pyabpoa

from svpg.util import sorted_nicely

class Candidate:
    def __init__(self, contig, start, end,  type, members, ref_seq='N', alt_seq='.', genotype='1/1', ref_reads=None, alt_reads=None, pan_known=None, detail_type=None, phase_list=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = type
        self.members = members
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.score = len(members)
        self.genotype = genotype
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads
        self.pan_known = pan_known
        self.detail_type = detail_type
        self.phase = phase_list

    def get_source(self):
        return (self.contig, self.start, self.end)

    def get_vcf_entry(self):
        contig, start, end = self.get_source()
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        if self.ref_reads != None and self.alt_reads != None:
            dp_string = str(self.ref_reads + self.alt_reads)
        else:
            dp_string = "."
        info_template = "SVTYPE={0};END={1};SVLEN={2};SUPPORT={3}"
        info_string = info_template.format(self.type, start if self.type == "INS" else end, start - end if self.type == "DEL" else end - start, self.score)
        if self.detail_type is not None:
            info_string += ";DETAILED_TYPE={0}".format(self.detail_type)
        if self.pan_known:
            info_string += ";PAN_KNOWN"
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
            chrom=contig,
            pos=start,
            id="PLACEHOLDERFORID",
            ref=self.ref_seq,
            alt=self.alt_seq,
            qual=self.score,
            filter="PASS" if len(filters) == 0 else ";".join(filters),
            info=info_string,
            format="GT:DP:AD",
            samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp=dp_string, ref=self.ref_reads if self.ref_reads != None else ".",
                                                   alt=self.alt_reads if self.alt_reads != None else "."))

class CandidateInversion(Candidate):
    def __init__(self, contig, start, source_direction, end, dest_direction, members,
                 support_fraction=".", genotype='1/1', ref_reads=None, alt_reads=None, detail_type=None):
        self.contig = contig
        #0-based source of the translocation (first base before the translocation)
        self.start = start
        self.direction = source_direction
        #0-based destination of the translocation (first base after the translocation)
        self.end = end
        self.dest_direction = dest_direction
        self.members = members
        self.score = len(members)
        self.type = "INV"
        self.support_fraction = support_fraction
        self.genotype = genotype
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads
        self.detail_type = detail_type

    def get_source(self):
        return (self.contig, self.start)

    def get_destination(self):
        return (self.contig, self.end)

    def get_vcf_entry(self):
        source_contig, source_start = self.get_source()
        dest_contig, dest_start = self.get_destination()

        if (self.direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "N[{contig}:{start}[".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "]{contig}:{start}]N".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig=dest_contig, start=dest_start)
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        info_template = "SVTYPE={0};SUPPORT={2}"
        info_string = info_template.format(self.type,  0, self.score)
        if self.detail_type is not None:
            info_string += ";DETAILED_TYPE={0}".format(self.detail_type)
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
            chrom=source_contig,
            pos=source_start,
            id="PLACEHOLDERFORID",
            ref="N",
            alt=alt_string,
            qual=self.score,
            filter="PASS" if len(filters) == 0 else ";".join(filters),
            info=info_string,
            format="GT:DP:AD",
            samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp='.', ref=".", alt="."))

class CandidateBreakend(Candidate):
    def __init__(self, source_contig, source_start, source_direction, dest_contig, dest_start, dest_direction, members,
                 support_fraction=".", genotype='1/1', ref_reads=None, alt_reads=None, detail_type=None):
        self.contig = source_contig
        #0-based source of the translocation (first base before the translocation)
        self.start = source_start
        self.direction = source_direction
        self.dest_contig = dest_contig
        #0-based destination of the translocation (first base after the translocation)
        self.dest_start = dest_start
        self.dest_direction = dest_direction
        self.members = members
        self.score = len(members)
        self.type = "BND"
        self.support_fraction = support_fraction
        self.genotype = genotype
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads
        self.detail_type = detail_type

    def get_source(self):
        return (self.contig, self.start)

    def get_destination(self):
        return (self.dest_contig, self.dest_start)

    def get_vcf_entry(self):
        source_contig, source_start = self.get_source()
        dest_contig, dest_start = self.get_destination()

        if (self.direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "N[{contig}:{start}[".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "]{contig}:{start}]N".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig=dest_contig, start=dest_start)
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        info_template = "SVTYPE={0};SUPPORT={2}"
        info_string = info_template.format(self.type,  0, self.score)
        if self.detail_type is not None:
            info_string += ";DETAILED_TYPE={0}".format(self.detail_type)
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
            chrom=source_contig,
            pos=source_start,
            id="PLACEHOLDERFORID",
            ref="N",
            alt=alt_string,
            qual=self.score,
            filter="PASS" if len(filters) == 0 else ";".join(filters),
            info=info_string,
            format="GT:DP:AD",
            samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp='.', ref=".", alt="."))

    def get_vcf_entry_reverse(self):
        source_contig, source_start = self.get_destination()
        dest_contig, dest_start = self.get_source()
        if (self.direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "N[{contig}:{start}[".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "]{contig}:{start}]N".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig=dest_contig, start=dest_start)
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        info_template = "SVTYPE={0};SUPPORT={1}"
        info_string = info_template.format(self.type, self.score)
        if self.detail_type is not None:
            info_string += ";DETAILED_TYPE={0}".format(self.detail_type)
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=source_contig,
                    pos=source_start,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt=alt_string,
                    qual=self.score,
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT:DP:AD",
                    samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp=".", ref=".", alt="."))

def _msa_consensus_for_cluster(cluster_seqs, aligner=None):
    if not cluster_seqs:
        return None
    if len(cluster_seqs) == 1:
        return cluster_seqs[0]

    if aligner is None:
        aligner = pyabpoa.msa_aligner()

    parts = sorted([(len(s), s, f'seq{i}') for i, s in enumerate(cluster_seqs)], reverse=True)
    _, seqs, names = zip(*parts)
    if not cluster_seqs:
        return []
    if aligner is None:
        aligner = pyabpoa.msa_aligner()

    aln_result = aligner.msa(list(seqs), out_msa=True, out_cons=True, max_n_cons=1)
    return aln_result.cons_seq

def consolidate_clusters_unilocal(clusters, ref_chrom, options, cons = False):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    min_sv_length, noseqs = options.min_sv_size, options.noseq
    max_sv_length = float('inf') if options.max_sv_size == -1 else options.max_sv_size
    ultra_sv_length = float('inf') if options.ultra_split_size == -1 else options.ultra_split_size
    aligner = pyabpoa.msa_aligner()
    repeat_pattern = re.compile(r'(A{20,}|T{20,}|(TC){20,}|(AG){20,})')

    consolidated_clusters = []
    for index, cluster in enumerate(clusters):
        svtype = cluster[0].type
        contig = cluster[0].get_source()[0]
        start = round(np.median([member.get_source()[1] for member in cluster]))

        members = [member.read_name for member in cluster]
        ref_seq = ref_chrom[max(start - 1, 0)]
        if svtype != "BND":
            end = round(np.median([member.svlen for member in cluster])) + start
            svlen = abs(end - start)
            if svlen > ultra_sv_length and len(members) < 10:
                continue
            if min_sv_length <= svlen <= max_sv_length:
                if not noseqs:
                    alt_seq = None
                    if svtype == "INS":
                        if not cons:
                            for member in cluster:
                                if member.svlen < svlen:
                                    continue
                                if member.alt_seq != "<INS>":
                                    alt_seq = ref_seq + member.alt_seq
                                else:
                                    alt_seq = "<INS>"
                                break
                        else:
                            seqs = []
                            for member in cluster:
                                if member.alt_seq != "<INS>":
                                    seqs.append(member.alt_seq)
                                else:
                                    alt_seq = "<INS>"
                                    break
                            if alt_seq != "<INS>" and svlen < 10000:
                                try:
                                    alt_seq = _msa_consensus_for_cluster(seqs, aligner=aligner)[0]
                                except Exception as e:
                                    alt_seq = seqs[0]
                            elif alt_seq != "<INS>":
                                alt_seq = seqs[0]
                    elif svtype == "DEL":
                        alt_seq = ref_seq
                        ref_seq = ref_chrom[max(start - 1, 0):end]
                        if options.read == "ont" and svlen < 100 and repeat_pattern.search(ref_seq):
                            continue
                    else:
                        ref_seq = "N"
                        alt_seq = f"<{svtype}>"
                else:
                    ref_seq = "N"
                    alt_seq = f"<{svtype}>"
                if svtype == "DUP":
                    consolidated_clusters.append(Candidate(contig, start, end,  svtype, members, ref_seq, alt_seq))
                elif svtype == "INV":
                    relative_direction = [inv.direction for inv in cluster]
                    common_relative_direction = Counter(relative_direction).most_common()
                    direction = common_relative_direction[0][0].split('_')[0]
                    detail_type = 'FOLDBACK_INV' if 'foldback' in cluster[0].direction else 'INV'
                    # consolidated_clusters.append(Candidate(contig, start, end, svtype, members, ref_seq, alt_seq, detail_type=detail_type))
                    if direction == 'left':
                        source_direction, dest_direction = 'fwd', 'rev'
                    else:
                        source_direction, dest_direction = 'rev', 'fwd'

                    consolidated_clusters.append(
                        CandidateBreakend(contig, start, source_direction, contig,
                                          end, dest_direction, members, detail_type=detail_type))

                else:  # INS,DEL
                    pan_known = True if cluster[0].node_ls else False
                    phase_list = [
                        getattr(member, 'phase')
                        for member in cluster
                        if hasattr(member, 'phase') and getattr(member, 'phase') is not None
                    ]
                    consolidated_clusters.append(Candidate(contig, start, end, svtype, members, ref_seq, alt_seq, pan_known=pan_known, phase_list=phase_list))
        else:
            dest_start = round(np.median([member.get_destination()[1] for member in cluster]))
            source_direction = max([member.source_direction for member in cluster],
                                   key=[member.source_direction for member in cluster].count)
            dest_direction = max([member.dest_direction for member in cluster],
                                 key=[member.dest_direction for member in cluster].count)
            # if contig == cluster[0].get_destination()[0]:
            #     detail_type = 'FOLDBACK_INV' if 'foldback' in cluster[0].source_direction else 'INV'
            #     consolidated_clusters.append(
            #         CandidateBreakend(contig, start, source_direction, contig, dest_start, dest_direction, members, detail_type=detail_type))
            #     consolidated_clusters.append(
            #         CandidateBreakend(contig, dest_start, source_direction, contig, start, dest_direction, members, detail_type=detail_type))
            consolidated_clusters.append(
                    CandidateBreakend(contig, start, source_direction, cluster[0].get_destination()[0], dest_start, dest_direction, members, detail_type='TRA'))

    return consolidated_clusters

def write_final_vcf(deletion_candidates,
                    novel_insertion_candidates,
                    duplication_candidates,
                    bnd_candidates,
                    contig_names,
                    contig_lengths,
                    options):
    types_to_output = [entry.strip() for entry in options.types.split(",")]
    vcf_output = open(os.path.join(options.working_dir, options.out), 'w')

    # Write header lines
    print("##fileformat=VCFv4.2", file=vcf_output)
    print("##fileDate={0}".format(time.strftime("%Y-%m-%d|%I:%M:%S%p|%Z|%z")), file=vcf_output)
    for contig_name, contig_length in zip(contig_names, contig_lengths):
        print("##contig=<ID={0},length={1}>".format(contig_name, contig_length), file=vcf_output)
    if "DEL" in types_to_output:
        print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    if "INS" in types_to_output:
        print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)

    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
          file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
          file=vcf_output)
    print("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting this variant\">",
          file=vcf_output)
    print("##INFO=<ID=DETAILED_TYPE,Number=1,Type=String,Description=\"Detailed type of the SV\">", file=vcf_output)
    if options.sub == 'call':
        print("##INFO=<ID=PAN_KNOWN,Number=0,Type=Flag,Description=\"Known structural variations in the pangenome\">",
              file=vcf_output)
    print("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">", file=vcf_output)
    print("##FILTER=<ID=hom_ref,Description=\"Genotype is homozygous reference\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">", file=vcf_output)
    print("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample", file=vcf_output)

    vcf_entries = []
    if "DEL" in types_to_output:
        for candidate in deletion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "DEL"))
    if "INS" in types_to_output:
        for candidate in novel_insertion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "INS"))
    if "DUP" in types_to_output:
        for candidate in duplication_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "DUP"))
    if "BND" in types_to_output or "INV" in types_to_output:
        pair_index = 0
        for candidate in bnd_candidates:
            entry1 = candidate.get_vcf_entry()
            entry2 = candidate.get_vcf_entry_reverse()

            # pair_id for linking two breakends
            pair_id = f"pair_{pair_index}"
            pair_index += 1

            vcf_entries.append(
                ((candidate.get_source()[0], candidate.get_source()[1], candidate.get_source()[1] + 1),
                 (entry1, pair_id), "BND"))
            vcf_entries.append(((candidate.get_destination()[0], candidate.get_destination()[1],
                                 candidate.get_destination()[1] + 1), (entry2, pair_id), "BND"))

    normalized = []
    for item in vcf_entries:
        loc = item[0]
        middle = item[1]
        svtype = item[2]

        if isinstance(middle, tuple) and len(middle) == 2:  # BND entry with pair_id
            entry_str, pair_id = middle
        else:   # non-BND
            entry_str, pair_id = middle, None

        normalized.append((loc, (entry_str, pair_id), svtype))

    entries_sorted = sorted_nicely(normalized)

    # assign variant IDs and handle MATEID for BNDs
    svtype_counter = defaultdict(int)
    pairid_to_svids = defaultdict(list)  # pair_id -> list of assigned svids
    output_buffer = []  # buffer for entries to write after MATEID assignment

    for loc, (entry_template, pair_id), svtype in entries_sorted:
        svtype_counter[svtype] += 1
        variant_id = f"SVPG.{svtype}.{svtype_counter[svtype]}"

        entry_assigned = entry_template.replace("PLACEHOLDERFORID", variant_id, 1)

        # record the svid for BND entries with pair_id
        if svtype == "BND" and pair_id is not None:
            pairid_to_svids[pair_id].append((variant_id, entry_assigned))
            # current svid and entry go to buffer for later MATEID filling
            output_buffer.append((pair_id, variant_id, entry_assigned, svtype))
        else:
            # non-BND entries go directly to output buffer
            output_buffer.append((None, variant_id, entry_assigned, svtype))

    #  if pair_id has two svids, fill in MATEID for both, then not fill MATEID
    for rec in output_buffer:
        pair_id, variant_id, entry_line, svtype = rec

        if svtype == "BND" and pair_id is not None:
            svid_list = pairid_to_svids.get(pair_id, [])
            if len(svid_list) == 2:
                svid_a, entry_a = svid_list[0]
                svid_b, entry_b = svid_list[1]

                mate_id = svid_b if variant_id == svid_a else svid_a

                # fill in MATEID
                if "MATE_PLACEHOLDER" in entry_line:
                    line_out = entry_line.replace("MATE_PLACEHOLDER", mate_id)
                else:
                    parts = entry_line.split("\t")
                    parts[7] = parts[7].rstrip() + f";MATEID={mate_id}"
                    line_out = "\t".join(parts)

                print(line_out, file=vcf_output)
            else:
                # only one breakend found for this pair_id, output without MATEID
                print(entry_line, file=vcf_output)
        else:
            print(entry_line, file=vcf_output)

    vcf_output.close()