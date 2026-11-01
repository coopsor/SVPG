import os
import subprocess
import re
from collections import Counter

import pysam.bcftools

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        total_written = 0
        currently_written = 0
        previous_label = ""
        nr_blocks = 0
        prev_n = False
        prev_break = True

        for line in infile:
            if line == "":
                continue
            if line.startswith(">"):
                previous_label = line.strip().split()[0]
                nr_blocks = 0
                prev_n = True
                continue
            for c in line:
                if c == '\n':
                    continue
                if c in ['n', 'N']:
                    prev_n = True
                    continue
                elif prev_n:
                    if not prev_break:
                        outfile.write("\n")
                    outfile.write(previous_label + "_" + str(nr_blocks) + "\n" + c)
                    total_written += 1
                    prev_break = False
                    currently_written = 1
                    nr_blocks += 1
                    prev_n = False

                else:
                    if currently_written % 70 == 0:
                        outfile.write("\n")
                        prev_break = True
                    outfile.write(c)
                    prev_break = False
                    currently_written += 1
                    total_written += 1
        if not prev_break:
            outfile.write("\n")
    print("Wrote " + str(total_written) + " non-N bases to output file.")


def filter_unsupported_variants(input_vcf, output_vcf):
    """
    Remove unsupported variant records before bcftools consensus.

    Currently filters:
      - Symbolic alleles, e.g. <INS>, <DEL>
      - Breakend alleles containing '[' or ']'
      - Spanning-deletion allele '*'
    """
    tmp_vcf = output_vcf + ".tmp.vcf"
    stats = Counter()

    result = subprocess.run(
            ["bcftools", "view", "-Ov", input_vcf],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
    
    with open(tmp_vcf, "w") as out:
        for line in result.stdout.splitlines():
            if line.startswith("#"):
                out.write(line + "\n")
                continue

            fields = line.split("\t")

            stats["total"] += 1
            alts = fields[4].split(",")
            unsupported = set()

            for alt in alts:
                if "[" in alt or "]" in alt:
                    unsupported.add("breakend")
                elif not re.fullmatch(r"[ACGT]+", alt, re.IGNORECASE):
                    unsupported.add("sequence_missing")

            if unsupported:
                stats["filtered"] += 1
                for variant_type in unsupported:
                    stats[variant_type] += 1
                continue

            out.write(line + "\n")

    # Convert to BGZF-compressed VCF, required for reliable tabix indexing.
    subprocess.run(
        [
            "bcftools", "view", "-Oz",
            "-o", output_vcf,
            tmp_vcf
        ],
        check=True
    )

    subprocess.run(
        ["bcftools", "index", "-f", output_vcf],
        check=True
    )
    print(
        ">>> VCF filtering summary: "
        f"total={stats['total']}, "
        f"filtered={stats['filtered']}, "
        f"sequence_missing={stats['sequence_missing']}, "
        f"breakend={stats['breakend']}"
    )


def augment_pipe(base_dir, ref_file, pan_file, output_file):
    os.chdir(base_dir)
    os.makedirs('tmp', exist_ok=True)

    print(">>> Merge VCF files using bcftools...")
    with open("filelist.tsv", "r") as f:
        input_vcfs = [line.strip() for line in f if line.strip()]

    if not input_vcfs:
        raise RuntimeError("No input VCF files found in filelist.tsv")
    if len(input_vcfs) == 1:
        subprocess.run(["bcftools", "view", "-Oz", "-o", "tmp/merged.vcf.gz", input_vcfs[0]], check=True)
    else:
        subprocess.run(["bcftools", "merge", "-l", "filelist.tsv", "--missing-to-ref", "--force-samples", "-Oz", "-o", "tmp/merged.vcf.gz"], check=True)
    pysam.tabix_index("tmp/merged.vcf.gz", preset="vcf", force=True)

    print(">>> Normalize the merged VCF file to ensure biallelic representation...")
    subprocess.run("bcftools norm -m- tmp/merged.vcf.gz  -Oz -o tmp/merged_biallelic.vcf.gz", shell=True, check=True)
    pysam.tabix_index("tmp/merged_biallelic.vcf.gz", preset="vcf", force=True)

    print(">>> Collapse the normalized VCF file using Truvari...")
    collapse_cmd = f"truvari collapse -i tmp/merged_biallelic.vcf.gz -o tmp/collapsed.tmp -c tmp/output.collapsed -f {ref_file}"
    subprocess.run(collapse_cmd, shell=True, check=True)
    subprocess.run("bcftools sort tmp/collapsed.tmp -W -Oz -o variants.vcf.gz", shell=True, check=True)

    print(">>> Generate the consensus sequence from the VCF file...")
    filter_unsupported_variants("variants.vcf.gz", "tmp/variants.filtered.vcf.gz")
    subprocess.run(f"bcftools consensus -f {ref_file} tmp/variants.filtered.vcf.gz | sed 's/^>/&augment_/' > tmp/cons.tmp", shell=True, check=True)

    process_fasta('tmp/cons.tmp', 'cons_noN.fa')
    subprocess.run("rm -rf tmp", shell=True, check=True)

    print(f">>> Construct a augment graph using Minigraph to {output_file}")
    subprocess.run(f"minigraph -cxggs -t128 {pan_file} cons_noN.fa > {output_file}", shell=True, check=True)