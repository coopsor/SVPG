import random
import vcf
import pysam

def read_variations(variation_file, ins_path, del_path):
    vcf_reader = vcf.Reader(open(variation_file, 'r'))
    vcf_writer = vcf.Writer(open(ins_path, 'w'), vcf_reader)
    with open(del_path, 'w') as f:
        for record in vcf_reader:
            values = [record.INFO.get('SVTYPE'), record.INFO.get('SVLEN')[0], record.INFO.get('END')]
            if values[0] == "DUP":
                record.INFO['SVTYPE'] = 'INS'
                bases = ['A', 'G', 'C', 'T']
                record.ALT = [''.join(random.choice(bases) for i in range(values[1]))]  # record.ALT is a list
                if values[2] != record.POS:
                    record.INFO['END'] = record.POS
                vcf_writer.write_record(record)
            elif values[0] == "DEL":
                # del_seq = ref.fetch(record.CHROM)[int(record.POS)-1:int(record.INFO['END'])]
                f.write(f'{record.CHROM}\t{record.POS - 1}\t{record.INFO["END"]}\n')
                # record.REF = ref.fetch(record.CHROM)[int(record.POS)-1:int(record.INFO['END'])]
                record.ALT = ['C']
                vcf_writer.write_record(record)
            # elif column[0] == "duplication":
            #     variations.append(variation.Insertion(values[2], values[1], values[0]))
            # elif column[0] == "translocation":
            #     variations.append(variation.Deletion(values[0], values[1]))
            #     variations.append(variation.Insertion(values[2], values[1], values[0]))
            else:
                pass


def del_ref(realsv, outsv, del_ref):
    vcf_reader = vcf.Reader(open(realsv, 'r'))
    vcf_writer = vcf.Writer(open(outsv, 'w'), vcf_reader)
    ref_genome = pysam.FastaFile('/data1/huheng/raresv/GRCh38_no_alt.fa')
    fasta_seq = {}
    with open(del_ref, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])

                region_seq = ref_genome.fetch(chrom)[start:end]
                key = f"{chrom}:{start}-{end}"
                fasta_seq[key] = region_seq

    for record in vcf_reader:
        if record.INFO.get('SVTYPE') == "DEL":
            record.REF = fasta_seq[''.join([record.CHROM, ':', str(record.POS - 1), '-', str(record.INFO['END'])])]
            # index += 1
        vcf_writer.write_record(record)


if __name__ == '__main__':
    variations = read_variations('/data1/huheng/raresv/filtered_sv.vcf', '/data1/huheng/raresv/realsv.vcf',
                                 '/data1/huheng/raresv/del_ref.bed')
    del_ref('/data1/huheng/raresv/realsv.vcf', '/data1/huheng/raresv/simsv.vcf', '/data1/huheng/raresv/del_ref.bed')
