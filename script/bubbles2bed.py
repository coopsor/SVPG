input_file = '/data/huheng/dataset/human/pangenome/GRCh38-90c.r518.bb.bed'
output_file = '/data/huheng/dataset/human/pangenome/GRCh38-90c.r518.sv.bed'
with open(input_file) as fin, open(output_file, 'w') as fout:
    for line in fin:
        fields = line.strip().split('\t')
        if len(fields) < 14:
            continue

        chrom = fields[0]
        start = int(fields[1])
        path = fields[11]
        for seq in fields[13:]:  # start from the 14th field
            end = start + len(seq)
            fout.write(f"{chrom}\t{start}\t{end}\t{path}\n")