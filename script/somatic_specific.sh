/home/huheng/SURVIVOR-master/Debug/SURVIVOR merge <(echo -e "/data1/huheng/HG008/variants.vcf\n/data1/huheng/HG008/variants_n.vcf") 1000 1 1 1 0 50 /data1/huheng/HG008/variants_merged.vcf

awk -F'\t' '
# If the line is a header (starts with #), print it as is
/^#/ {
    print $0;
    next;
}

{
    # Extract genotype fields from the last two columns
    split($(NF-1), format1, ":");  # Second-to-last column (e.g., sample1)
    split($NF, format2, ":");      # Last column (e.g., sample2)

    gt1 = format1[1];  # Genotype of sample1
    gt2 = format2[1];  # Genotype of sample2

    # Print variants where sample1 is not homozygous reference or missing,
    # and sample2 is homozygous reference or missing
    if (gt1 != "0/0" && gt1 != "./." && (gt2 == "0/0" || gt2 == "./.")) {
        print $0;
    }
}
' /data1/huheng/HG008/variants_merged.vcf > /data1/huheng/HG008/variants_specific.vcf

