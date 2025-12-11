#!/usr/bin/env python3
import sys
import re
from collections import defaultdict
import os

"""
Evaluate SV calls by type (INDEL, INV, TRA) between a base VCF and one or more comparison VCFs.
"""

def bedtools_filter(vcf_path, bed_path):
    """using bedtools intersect filter VCF，only retain SV within bed regions"""
    tmp_vcf = os.path.abspath(os.path.basename(vcf_path).replace('.vcf', '_filtered.vcf'))
    bed_regions = {}
    with open(bed_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.strip().split()[:3]
            start, end = int(start), int(end)
            if chrom not in bed_regions:
                bed_regions[chrom] = []
            bed_regions[chrom].append((start, end))
    for chrom in bed_regions:
        bed_regions[chrom].sort()

    with open(vcf_path) as vin, open(tmp_vcf, "w") as vout:
        for line in vin:
            if line.startswith("#"):
                vout.write(line)
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]

            end = None
            svlen = None
            for kv in info.split(";"):
                if kv.startswith("END="):
                    try:
                        end = int(kv.split("=")[1])
                    except ValueError:
                        pass
                elif kv.startswith("SVLEN="):
                    try:
                        svlen = abs(int(kv.split("=")[1]))
                    except ValueError:
                        pass

            if end is None and svlen is not None:
                end = pos + svlen
            if end is None:
                end = pos

            if chrom not in bed_regions:
                continue

            in_bed = False
            for start_b, end_b in bed_regions[chrom]:
                if end_b < pos:
                    continue
                if start_b > end:
                    break
                if not (end < start_b or pos > end_b):
                    in_bed = True
                    break

            if in_bed:
                vout.write(line)
    return tmp_vcf


def load_bed(bed_path):
    bed_regions = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            bed_regions[chrom].append((start, end))
    for chrom in bed_regions:
        bed_regions[chrom].sort()
    return bed_regions


def in_bed(chrom, start, end, bed_regions):
    if chrom not in bed_regions:
        return False

    for bed_start, bed_end in bed_regions[chrom]:
        if end < bed_start:
            break
        if bed_start <= start <= bed_end or bed_start <= end <= bed_end:
            return True

    return False

def chr_to_sort_key(chr_name):
    if chr_name.startswith("chr"):
        chr_name = chr_name[3:]
    if chr_name.isdigit():
        return int(chr_name)
    elif chr_name in ["X", "Y", "M", "MT"]:
        return {"X": 23, "Y": 24, "M": 25, "MT": 25}[chr_name]
    else:
        return None

def parse_vcf(vcf_path, bed_regions=None):
    sv_records = defaultdict(list)

    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, id, ref, alt, qual, filter, info = fields[:8]
            if filter != 'PASS':
                continue
            pos = int(pos)
            svlen = None
            line_upper = line.upper()

            # --- find SV type ---
            svtype = None
            no_bnd_type = ['DEL', 'INS', 'DUP', 'DUP:TANDEM', 'INV', 'FOLDBACK', 'TRA']
            matches = re.findall(
                r'\b(SVTYPE|EVENTTYPE|DETAILED_TYPE|SV_CLASS|SV_TYPE|TYPE)=([\w:]+)',
                info, re.IGNORECASE
            )
            if matches:
                raw_types = [m[1].upper() for m in matches if not m[1].isdigit()]
                for candidate in raw_types:
                    for t in no_bnd_type:
                        if t in candidate:
                            svtype = t
                            break
            if not svtype:
                for t in no_bnd_type:
                    if re.search(rf'\b{t}\b', line_upper, re.IGNORECASE):
                        svtype = t
                        break

            svtype_map = {
                'INS': 'INDEL',
                'DEL': 'INDEL',
                'DUP': 'INDEL',
                'DUP:TANDEM': 'INDEL',
                'INV': 'INV',
                'FOLDBACK': 'INV',
                'TRA': 'TRA'
            }
            if svtype:
                svtype = svtype_map.get(svtype, svtype)

            end = None
            if 'END=' in info:
                try:
                    end = re.search(r'(?i)\bEND=([^;]+)', info).group(1)
                except:
                    pass
            if 'SVLEN=' in info:
                try:
                    svlen = abs(int(re.search(r'SVLEN=(-?\d+)', info).group(1)))
                except:
                    pass

            # --- BND handling ---
            if ('[' in alt or ']' in alt) or ('[' in info or ']' in info):
                pattern = []

                if isinstance(alt, str):
                    pattern = re.findall(r'[\[\]](chr[\w]+):(\d+)[\[\]]', alt)

                # END=N[chr2:8888]
                if not pattern and isinstance(end, str):
                    pattern = re.findall(r'[\[\]](chr[\w]+):(\d+)[\[\]]', end)

                if pattern:
                    chr2, pos2 = pattern[0]
                    pos2 = int(pos2)
                    if chr_to_sort_key(chrom) is None or chr_to_sort_key(chr2) is None:
                        continue

                    if chrom != chr2:
                        sv_records['TRA'].append({
                            'chrom': chrom,
                            'chrom_next': chr2,
                            'start': pos,
                            'end': pos2
                        })
                    else:
                        if svtype == 'INV':
                            sv_records['INV'].append({
                                'chrom': chrom,
                                'chrom_next': chr2,
                                'start': pos,
                                'end': pos2
                            })

                    continue

            if not svtype or svtype == 'INV':
                if end:
                    sv_records['INV'].append({
                        'chrom': chrom,
                        'chrom_next': chrom,
                        'start': pos,
                        'end': int(end)
                    })
                continue

            if not svlen:
                if end:
                    svlen = abs(int(end) - pos)
                else:
                    continue

            # --- non-BND ---
            if svlen < 50:
                continue
            # if not bed_regions or in_bed(chrom, pos, end, bed_regions):
            sv_records[svtype].append({
                'chrom': chrom,
                'start': pos,
                'len': svlen
            })

    return sv_records


def match_sv(base_list, comp_list, svtype):
    matched = set()
    unmatched = []
    false_positives = []
    tp = 0

    for i, b in enumerate(base_list):
        found = False
        for j, c in enumerate(comp_list):

            if svtype == 'TRA' or svtype == 'INV':
                if b['chrom'] == c['chrom'] and b['chrom_next'] == c['chrom_next'] and abs(b['start'] - c['start']) <= 1000 and abs(b['end'] - c['end']) <= 1000:
                    key_contig, key_start = b['chrom'], b['start']
                    tp += 1
                    matched.add(f"{key_contig}-{key_start}")
                    found = True
                    break
            else:
                if c['chrom'] != b['chrom']:
                    continue
                key_contig, key_start = b['chrom'], b['start']
                if abs(b['start'] - c['start']) <= 500:
                    try:
                        ratio = min(b['len'], c['len']) / max(b['len'], c['len'])
                    except:
                        ratio = 0
                    if ratio >= 0.7:
                        tp += 1
                        matched.add(f"{key_contig}-{key_start}")
                        found = True
                        break

        if not found:
            unmatched.append(b)


    matched_keys = {m.split('-')[0] + ':' + m.split('-')[1] for m in matched}
    for c in comp_list:
        comp_key = f"{c['chrom']}:{c['start']}"
        if not any(abs(c['start'] - b['start']) <= 500 and c['chrom'] == b['chrom'] for b in base_list):
            false_positives.append(c)

    fp = len(false_positives)
    fn = len(unmatched)

    # print matching summary
    # print(f"\n=== {svtype} 匹配结果 ===")
    # print(f"TP: {tp}, FP: {fp}, FN: {fn}")
    # if unmatched:
    #     print(f"\n[FN]: ({len(unmatched)} 个)")
    #     for b in unmatched[:10]:  # 只打印前10个
    #         print(f"  {b['chrom']}:{b['start']}-{b.get('end', '?')}  len={b.get('len', '?')}")
    #     if len(unmatched) > 10:
    #         print(f"  ...  {len(unmatched) - 10} not shown")
    #
    # if false_positives:
    #     print(f"\n[FP]: ({len(false_positives)} 个)")
    #     for c in false_positives[:10]:
    #         print(f"  {c['chrom']}:{c['start']}-{c.get('end', '?')}  len={c.get('len', '?')}")
    #     if len(false_positives) > 10:
    #         print(f"  ... {len(false_positives) - 10} not shown")

    return tp, fp, fn, matched, (unmatched, false_positives)


def calc_metrics(tp, fp, fn):
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
    return prec, rec, f1


def main():
    if len(sys.argv) < 4:
        print("Usage: python eval_svtype_split.py base.vcf comp.vcf")
        sys.exit(1)
    else:
        print(' '.join(sys.argv))
    base_vcf = sys.argv[1]

    last_arg = sys.argv[-1]
    has_bed = last_arg.endswith(".bed")

    if has_bed:
        bed_path = last_arg
        comp_vcfs = sys.argv[2:-1]
    else:
        bed_path = None
        comp_vcfs = sys.argv[2:]

    if bed_path:
        print(f"[INFO] Filtering VCFs by BED using bedtools...")
        base_vcf = bedtools_filter(base_vcf, bed_path)
        comp_vcfs = [bedtools_filter(vcf, bed_path) for vcf in comp_vcfs]

    base_sv = parse_vcf(base_vcf)

    svtypes = ['INDEL', 'INV', 'TRA']
    # {svtype: {tool: (tp, fp, fn, prec, rec, f1)}}
    summary = {svt: {} for svt in svtypes}
    summary["ALL"] = {}
    for index, comp_vcf in enumerate(comp_vcfs):
        print(f"\n=== Evaluating: {comp_vcf} ===")
        tool_index = f"Tool_{index+1}"
        comp_sv = parse_vcf(comp_vcf)
        total_tp = total_fp = total_fn = 0

        print("\t".join(['SVTYPE', 'TP', 'FP', 'FN', 'Precision', 'Recall', 'F1']))
        for svt in svtypes:
            tp, fp, fn, matched, unmatched = match_sv(base_sv.get(svt, []), comp_sv.get(svt, []), svt)
            prec, rec, f1 = calc_metrics(tp, fp, fn)
            total_tp += tp
            total_fp += fp
            total_fn += fn
            print(f"{svt}\t{tp}\t{fp}\t{fn}\t{prec:.4f}\t{rec:.4f}\t{f1:.4f}")
            summary[svt][tool_index] = (tp, fp, fn, prec, rec, f1)
            # print(unmatched)

        prec, rec, f1 = calc_metrics(total_tp, total_fp, total_fn)
        print(f"ALL\t{total_tp}\t{total_fp}\t{total_fn}\t{prec:.4f}\t{rec:.4f}\t{f1:.4f}")
        summary["ALL"][tool_index] = (total_tp, total_fp, total_fn, prec, rec, f1)

    print("\n=== Summary across all SV type ===")
    header = ["Tool", "TP", "FP", "FN", "Precision", "Recall", "F1"]
    print("\t".join(header))
    for svt in summary:
        print(f"Summary of {svt}")
        for tool, vals in summary[svt].items():
            tp, fp, fn, prec, rec, f1 = vals
            print(f"{tool}\t{tp}\t{fp}\t{fn}\t{prec:.4f}\t{rec:.4f}\t{f1:.4f}")

if __name__ == "__main__":
    main()
