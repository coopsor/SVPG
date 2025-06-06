import pandas as pd

df = pd.read_excel('/home/huheng/code/SVdetect/SVScan/13059_2022_2816_MOESM4_ESM.xlsx')


vcf_header = """##fileformat=VCFv4.2
##source=Excel2VCF
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

vcf_lines = []
for _, row in df.iterrows():
    try:
        chrom = row['Chrom1']
        # if not row['Pos1']:
        #     break
        pos = int(row['Pos1'])
        sv_id = row['SV_id']
        sv_type = row['SV_type']
        end = int(row['Pos2'])
        if sv_type == 'DEL':
            svlen = int(row['SV_Size or breakpoints distance'])*(-1)
        elif sv_type == 'INS':
            svlen = int(row['SV_Size or breakpoints distance'])
        else:
            continue

        alt = f"<{sv_type}>"
        info = f"SVTYPE={sv_type};END={end};SVLEN={svlen}"

        format_field = "GT:DP:AD"
        sample_field = "1/1:.:.,."

        line = f"{chrom}\t{pos}\t{sv_id}\tN\t{alt}\t.\tPASS\t{info}\t{format_field}\t{sample_field}"
        vcf_lines.append(line)
    except Exception as e:
        break

with open("/data/huheng/dataset/human/Tumor/output4.vcf", "w") as f:
    f.write(vcf_header)
    f.write("\n".join(vcf_lines))
