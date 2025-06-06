import numpy as np
import time
import os.path
from collections import defaultdict

from svpg.util import sorted_nicely

class Candidate:
    def __init__(self, contig, start, end, sv_from, type, members, ref_seq='N', alt_seq='.', genotype='1/1', ref_reads=None, alt_reads=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.sv_from = sv_from
        self.type = type
        self.members = members
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.score = len(members)
        self.genotype = genotype
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads

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
        info_string = info_template.format(self.type, end if self.type == "DEL" else start, start - end if self.type == "DEL" else end - start, self.score)
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

def consolidate_clusters_unilocal(clusters, ref_chrom, options, fasta_file=None):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    min_sv_length, noseqs = options.min_sv_size, options.noseq
    max_sv_length = float('inf') if options.max_sv_size == -1 else options.max_sv_size

    consolidated_clusters = []
    for index, cluster in enumerate(clusters):
        svtype = cluster[0].type
        sv_from = [member.signature for member in cluster]
        members = [member.read_name for member in cluster]
        start = round(np.median([member.get_source()[1] for member in cluster]))
        end = round(np.median([member.svlen for member in cluster])) + start
        svlen = abs(end - start)

        if min_sv_length <= svlen <= max_sv_length:
            if not noseqs:
                try:
                    ref_base = ref_chrom[max(start - 1, 0)]
                    ref_seq = ref_base if svtype == "INS" else ref_chrom[max(start - 1, 0):end]
                except IndexError:
                    continue

                if svtype == "INS":
                    for member in cluster:
                        if member.svlen < svlen:
                            continue
                        if member.alt_seq != "<INS>":
                            alt_seq = ref_base + member.alt_seq
                        else:
                            if fasta_file:
                                alt_seq = ref_base
                            else:
                                alt_seq = "<INS>"
                        break
                else:
                    alt_seq = ref_base if ref_chrom else "<DEL>"
            else:
                ref_seq = "N"
                alt_seq = "<INS>" if svtype == "INS" else "<DEL>"

            consolidated_clusters.append(Candidate(cluster[0].get_source()[0], start, end, sv_from, svtype, members, ref_seq, alt_seq))

    return consolidated_clusters

def write_final_vcf(deletion_candidates,
                    novel_insertion_candidates,
                    contig_names,
                    contig_lengths,
                    types_to_output,
                    working_dir,
                    outfile):
    vcf_output = open(os.path.join(working_dir, outfile), 'w')

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
    print("##FILTER=<ID=hom_ref,Description=\"Genotype is homozygous reference\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">", file=vcf_output)
    print("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample", file=vcf_output)

    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
    if "DEL" in types_to_output:
        for candidate in deletion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "DEL"))
    if "INS" in types_to_output:
        for candidate in novel_insertion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "INS"))

    # Sort and write entries to VCF
    svtype_counter = defaultdict(int)
    for source, entry, svtype in sorted_nicely(vcf_entries):
        variant_id = "SVPG.{svtype}.{number}".format(svtype=svtype, number=svtype_counter[svtype] + 1)
        entry_with_id = entry.replace("PLACEHOLDERFORID", variant_id, 1)
        svtype_counter[svtype] += 1
        print(entry_with_id, file=vcf_output)

    vcf_output.close()
