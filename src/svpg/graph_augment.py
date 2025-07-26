import os
import subprocess

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


def augment_pipe(base_dir, ref_file, pan_file, output_file):
    os.chdir(base_dir)
    os.makedirs('tmp', exist_ok=True)

    print(">>> Merge VCF files using bcftools...")
    pysam.bcftools.merge("-l", "filelist.tsv", "--missing-to-ref", "--force-samples", "-Oz", "-o", "tmp/merged.vcf.gz", force=True, catch_stdout=False)
    pysam.tabix_index("tmp/merged.vcf.gz", preset="vcf", force=True)

    print(">>> Normalize the merged VCF file to ensure biallelic representation...")
    subprocess.run("bcftools norm -m- tmp/merged.vcf.gz  -Oz -o tmp/merged_biallelic.vcf.gz", shell=True, check=True)
    pysam.tabix_index("tmp/merged_biallelic.vcf.gz", preset="vcf", force=True)

    print(">>> Collapse the normalized VCF file using Truvari...")
    collapse_cmd = f"truvari collapse -i tmp/merged_biallelic.vcf.gz -o tmp/collapsed.tmp -c tmp/output.collapsed -f {ref_file}"
    subprocess.run(collapse_cmd, shell=True, check=True)
    subprocess.run("bcftools sort tmp/collapsed.tmp -W -Oz -o variants.vcf.gz", shell=True, check=True)

    print(">>> Generate the consensus sequence from the VCF file...")
    subprocess.run(f"bcftools consensus -f {ref_file} variants.vcf.gz | sed 's/^>/&augment_/' > tmp/cons.tmp", shell=True, check=True)

    process_fasta('tmp/cons.tmp', 'cons_noN.fa')
    subprocess.run("rm -rf tmp", shell=True, check=True)

    print(f">>> Construct a augment graph using Minigraph to {output_file}")
    subprocess.run(f"minigraph -cxggs -t128 {pan_file} cons_noN.fa > {output_file}", shell=True, check=True)