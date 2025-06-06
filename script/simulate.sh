# Extract header lines from the original VCF
grep '^#' raw_raresv.vcf > filtered_sv.vcf

# Filter for high-quality rare DEL/DUP SVs found in only one sample or one family
awk -F'\t' '$7 == "PASS" && ($8 ~ /SVTYPE=DEL/ || $8 ~ /SVTYPE=DUP/) && ($8 ~ /NSAMP=1/ || $8 ~ /NFAM=1/)' raw_raresv.vcf >> filtered_sv.vcf

# Convert VCF to BED format for overlap testing with pangenome bubbles,
# representing the rarity level of structural variants
awk -F'\t' 'BEGIN {OFS="\t"} !/^#/ {split($8, info, ";"); start=$2; end=""; svtype=""; for (i in info) {if (info[i] ~ /^END=/) end=substr(info[i], 5); if (info[i] ~ /^SVTYPE=/) svtype=substr(info[i], 8)}; print $1, start, end, svtype}' filtered_sv.vcf > output.bed

# Calculate overlap (Jaccard index) between filtered SVs and pangenome graph bubbles
bedtools jaccard -a output.bed -b /data/huheng/pangenome/GRCh38-90c.r518.bb.bed

# Add genotype (GT) fields to VCF for downstream tools like truvari
awk 'BEGIN {OFS="\t"}
  /^#/ {
      if ($0 ~ /^#CHROM/) {
          print $0, "FORMAT", "SAMPLE"
      } else {
          print $0
      }
  }
  !/^#/ {
      print $0, "GT", "1/1"
  }' /data1/huheng/pangenome_graph/filtered_sv.vcf > /data1/huheng/pangenome_graph/filtered_gt.vcf

# Generate simulated SV VCF using simulated_Add_alt, then construct consensus sequence with bcftools
bcftools consensus -f GRCh38_no_alt.fa simsv.vcf.gz -o simsv.fa 2 > skip.log

# Filter skipped positions reported by bcftools to construct a clean SV benchmark set
grep "overlaps with another variant, skipping" skip.log | awk '{print $3}' > positions.txt

awk 'BEGIN {
    while ((getline < "positions.txt") > 0) {
        split($0, pos, ":")
        positions[pos[1] ":" pos[2]] = 1
    }
    close("positions.txt")
}
{
    if (!($1 ":" $2 in positions)) print $0
}' filtered_gt.vcf > new_filtered.vcf

# Simulate sequencing reads from the consensus genome using PBSIM
pbsim --strategy wgs --method qshmm --qshmm R103.model --depth 20 --genome simsv.fa
