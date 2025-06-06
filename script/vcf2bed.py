import vcf

variation_file = '/home/huheng/code/PanSV/eval_raresv_svpg/tp.vcf'
vcf_reader = vcf.Reader(open(variation_file, 'r'))
bed_file = open('/home/huheng/code/PanSV/eval_raresv_svpg/tp.bed', 'w')
for record in vcf_reader:
    if record.INFO.get('SVTYPE') not in ['INS', 'DEL']:
        continue
    values = [record.INFO.get('SVTYPE'), record.INFO.get('SVLEN'), record.INFO.get('END')]
    if values[0] == 'INS':
        end = record.POS + values[1][0]  # values[1] is a list, take the first element
    else:
        end = values[2]
    bed_file.write(f'{record.CHROM}\t{record.POS}\t{end}\t{values[0]}\n')


