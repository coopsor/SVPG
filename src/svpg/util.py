import os
import re

class gfaNode:
    def __init__(self, name="", sequence="", length=0, contig="", offset=0, sr=0):
        self.name = name
        self.sequence = sequence
        self.len = length
        self.contig = contig
        self.offset = offset
        self.sr = sr

def parse_contig_from_sn(tag: str) -> str:
    """
    Parse contig name from GFA SN:Z: tag.
    Supports:
      - SN:Z:CHM13#0#chr1
      - SN:Z:id=CHM13v2|chr1
      - SN:Z:chr1
    """
    assert tag.startswith("SN:Z:"), f"Unexpected SN tag: {tag}"

    value = tag[5:]  # remove 'SN:Z:'

    # case 1: contains '|', take last field
    if '|' in value:
        return value.split('|')[-1]

    # case 2: contains '#', take last field
    if '#' in value:
        return value.split('#')[-1]

    # case 3: exit
    return value


def read_gfa(ref_graph):
    """Parse a GFA file to extract node and edge information."""
    gfa_node, gfa_edge = {}, {}
    with open(ref_graph, "r") as fp:
        for line in fp:
            tokens = line.strip().split('\t')
            if tokens[0] == "S":  # Sequence line
                name = tokens[1]
                sequence = tokens[2]
                try:
                    length = int(tokens[3][5:])  # e.g., LN:i:1234 -> 1234
                    contig = parse_contig_from_sn(tokens[4])  # e.g., SN:Z:chr1 -> chr1
                    offset = int(tokens[5][5:])
                    sr = int(tokens[6][5:])  # e.g., SR:i:1
                except (IndexError, ValueError):
                    raise ValueError(f"Invalid S-line format in GFA: {line}")

                gfa_node[name] = gfaNode(name, sequence, length, contig, offset, sr)

            # elif tokens[0] == "L":  # Link/edge line
            #     gfa_edge.setdefault(tokens[1], []).append(line)

    return gfa_node

def analyze_cigar_indel(tuples, min_length, is_gaf=False):
    """
    Parses CIGAR tuples and returns indels with length >= min_length.
    If is_gaf=True, expects GAF-style tuples (length, op).
    Otherwise, expects BAM-style tuples (op, length).
    """
    pos_ref = 0
    pos_read = 0  #
    indels = []

    # BAM CIGAR code to string mapping
    bam_op_map = {
        0: 'M',  # match
        1: 'I',  # insertion
        2: 'D',  # deletion
        4: 'S',  # soft clip
        7: '=',  # match
        8: 'X',  # mismatch
    }

    for t in tuples:
        if is_gaf:
            length, op = t
        else:
            op_code, length = t
            op = bam_op_map.get(op_code, None)
            if op is None:
                continue  # skip unknown op codes

        if op == 'M' or op == '=' or op == 'X':
            pos_ref += length
            pos_read += length
        elif op == 'I':
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS"))
            pos_read += length
        elif op == 'D':
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL"))
            pos_ref += length
        elif op == 'S':
            pos_read += length

    return indels


def merge_cigar(sigs, max_merge=500):
    """
       Merge nearby SV signatures of the same type if they are within max_merge distance
       and each is longer than min_indel_length.
       """
    i = 0
    while i < len(sigs) - 1:
        diff = abs(sigs[i + 1].start - sigs[i].start) if sigs[i].type == 'INS' else abs(sigs[i + 1].start - sigs[i].end)
        if (sigs[i].type == sigs[i + 1].type and diff <= max_merge):
            sigs[i].end += sigs[i + 1].svlen
            sigs[i].svlen += sigs[i + 1].svlen
            sigs.pop(i + 1)
        else:
            i += 1

    return sigs

def chr_to_sort_key(chr_name):
    if chr_name.startswith("chr"):
        chr_name = chr_name[3:]
    if chr_name.isdigit():
        return int(chr_name)
    elif chr_name in ["X", "Y", "M", "MT"]:
        return {"X": 23, "Y": 24, "M": 25, "MT": 25}[chr_name]
    else:
        return None

def sorted_nicely(vcf_entries):
    """ Sort the given vcf entries (in the form ((contig, start, end), vcf_string, sv_type)) in the way that humans expect.
        e.g. chr10 comes after chr2
        Algorithm adapted from https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    tuple_key = lambda entry: (alphanum_key(str(entry[0][0])), entry[0][1], entry[0][2])
    return sorted(vcf_entries, key=tuple_key)


def find_sequence_file(entry):
    """Find a sequence file in the given entry's path with common suffixes."""
    common_suffixes = ['.fa', '.fasta', '.fna', '.fastq', '.fq',
                       '.fa.gz', '.fasta.gz', '.fna.gz', '.fastq.gz', '.fq.gz']

    for suffix in common_suffixes:
        filename = f"{entry.name}{suffix}"
        full_path = os.path.join(entry.path, filename)
        if os.path.exists(full_path):
            return suffix

    return None
