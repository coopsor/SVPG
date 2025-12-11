#!/usr/bin/env python3
import sys

'''
Get cancer-specific SVs from vcf1(tumor.vcf) compared to vcf2(normal.vcf)
'''

def get_svtype(line):
    line_upper = line.upper()
    if 'SVTYPE=INS' in line_upper or 'EVENTTYPE=INS' in line_upper or 'DETAILED_TYPE=INS' in line_upper:
        return 'INS'
    elif 'SVTYPE=DEL' in line_upper or 'EVENTTYPE=DEL' in line_upper or 'DETAILED_TYPE=DEL' in line_upper:
        return 'DEL'
    elif 'SVTYPE=DUP' in line_upper or 'EVENTTYPE=DUP' in line_upper or 'DETAILED_TYPE=DUP' in line_upper or 'DUP' in line_upper:
        return 'DUP'
    elif 'SVTYPE=INV' in line_upper or 'EVENTTYPE=INV' in line_upper or 'DETAILED_TYPE=INV' in line_upper or 'INV' in line_upper or 'FOLDBACK' in line_upper:
        return 'INV'
    else:
        return 'TRA'

def parse_vcf(vcf_file):
    sv_records = {}
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]
            # m = re.search(r"SVTYPE=([^;]+)", info)
            # if not m:
            #     continue
            # svtype = m.group(1)
            line_upper = line.upper()
            # --- 判断SV类型 ---
            svtype = get_svtype(line)
            sv_records.setdefault(chrom, {}).setdefault(svtype, []).append(line)
    return sv_records

def is_nearby(pos, line_list, threshold=1000):
    pos_list = [int(rec.split('\t')[1])for rec in line_list]
    for p, line in zip(pos_list, line_list):
        if abs(p - pos) <= threshold:
            return line
    return False


def merge_vcf_records(vcf1_file, vcf2_file, output_file):
    sv2 = parse_vcf(vcf2_file)
    with open(vcf1_file) as f1, open(output_file, 'w') as out:
        for line in f1:
            if line.startswith('#'):
                out.write(line)
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            svtype = get_svtype(line)

            if chrom in sv2 and svtype in sv2[chrom]:
                line_other = is_nearby(pos, sv2[chrom][svtype], 500)
                if line_other:
                    out.write(line_other)
                    continue


def merge_vcf_specific_records(vcf1, vcf2, output):
    sv2 = parse_vcf(vcf2)

    with open(vcf1) as f1, open(output, 'w') as out:
        for line in f1:
            if line.startswith('#'):
                out.write(line)
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]
            # m = re.search(r"SVTYPE=([^;]+)", info)
            # if not m:
            #     continue
            # svtype = m.group(1)

            svtype = get_svtype(line)

            if chrom in sv2 and svtype in sv2[chrom]:
                if is_nearby(pos, sv2[chrom][svtype], 1000):
                    continue

            out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_specific.py <vcf1> <vcf2> <output>")
        sys.exit(1)

    vcf1, vcf2, output = sys.argv[1], sys.argv[2], sys.argv[3]
    merge_vcf_specific_records(vcf1, vcf2, output)