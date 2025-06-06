#!/bin/bash

vcfa="/home/huheng/code/PanSV/ont_20x/var_.vcf"
vcfb="/home/huheng/code/PanSV/ont_20x_/var_.vcf"

echo -e "$vcfa\n$vcfb" > /home/huheng/SURVIVOR-master/Debug/sample_files

# intersect
/home/huheng/SURVIVOR-master/Debug/SURVIVOR merge /home/huheng/SURVIVOR-master/Debug/sample_files 1000 2 1 1 0 50 /home/huheng/SURVIVOR-master/Debug/sample_merged.vcf

# union
/home/huheng/SURVIVOR-master/Debug/SURVIVOR merge /home/huheng/SURVIVOR-master/Debug/sample_files 1000 1 1 1 0 50 /home/huheng/SURVIVOR-master/Debug/sample_merged_.vcf

# count non-header lines in both VCF files
count1=$(grep -vc "^#" /home/huheng/SURVIVOR-master/Debug/sample_merged.vcf)
count2=$(grep -vc "^#" /home/huheng/SURVIVOR-master/Debug/sample_merged_.vcf)
echo "count1: $count1"
echo "count2: $count2"

if [ "$count2" -ne 0 ]; then
    ratio=$(echo "scale=3; $count1 / $count2" | bc)
    printf "Ratio of non-header lines: %0.3f\n" "$ratio"
else
    echo "The second VCF file has no non-header lines, division by zero."
fi
