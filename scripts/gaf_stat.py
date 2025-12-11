import re
import sys
import numpy as np
from collections import defaultdict, Counter

""""
Analyze GAF alignment files and generate statistics.
"""

class Gaf:
    def __init__(self):
        self.query_name = ""
        self.query_length = 0
        self.query_start = 0
        self.query_end = 0
        self.strand = ""
        self.path = ""
        self.path_length = 0
        self.path_start = 0
        self.path_end = 0
        self.mapping_quality = 0
        self.is_primary = True
        self.cigar = ""
        self.ds = ""


def parse_gaf_line(tokens, gfa_node):
    """Parse a single GAF line"""
    gafline = Gaf()

    gafline.query_name = tokens[0]
    gafline.query_length = int(tokens[1])
    gafline.query_start = int(tokens[2])
    gafline.query_end = int(tokens[3])

    gafline.path = re.findall(r'[<>][^<>]+', tokens[5])
    try:
        gafline.strand = '+' if gafline.path[0][0] == '>' else '-'
    except:
        pass

    # Ensure gfa_node dictionary contains this path node
    try:
        gafline.contig = gfa_node[gafline.path[0][1:]].contig
        gafline.offset = gfa_node[gafline.path[0][1:]].offset
    except:
        # Set default values if contig and offset are unavailable
        gafline.contig = "unknown"
        gafline.offset = 0

    gafline.path_length = int(tokens[6])
    gafline.path_start = int(tokens[7])
    gafline.path_end = int(tokens[8])
    gafline.mapping_quality = int(tokens[11])
    gafline.is_primary = True

    for tok in tokens:
        if "tp:A:" in tok:
            if tok[5:7] != "P":
                gafline.is_primary = False
        if "cg:Z:" in tok:
            gafline.cigar = tok[5:]
        if "ds:Z:" in tok:
            gafline.ds = tok[6:]

    return gafline


def analyze_gaf_file(gaf_file, gfa_node_dict):
    """
    Analyze a GAF file and generate detailed statistics.

    Args:
        gaf_file (str): Path to the GAF file.
        gfa_node_dict (dict): Dictionary of GFA node information.

    Returns:
        dict: A dictionary of various statistics.
    """
    stats = {
        'total_alignments': 0,
        'primary_alignments': 0,
        'secondary_alignments': 0,
        'unique_reads': set(),
        'mapped_bases': 0,
        'total_query_bases': 0,
        'mapq_sum': 0,
        'mapq_by_range': {
            '0': 0,
            '1-10': 0,
            '11-20': 0,
            '21-30': 0,
            '31-40': 0,
            '41-50': 0,
            '51-60': 0,
            '>60': 0
        },
        'reads_with_multiple_alignments': 0,
        'alignment_lengths': [],
        'query_coverage': [],
        'multi_path_alignments': 0,
        'strand_counts': {'forward': 0, 'reverse': 0}
    }

    read_alignments = defaultdict(list)

    try:
        with open(gaf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip comment lines

                tokens = line.strip().split('\t')
                if len(tokens) < 12:
                    continue  # Skip malformed lines

                gaf = parse_gaf_line(tokens, gfa_node_dict)
                stats['total_alignments'] += 1
                stats['unique_reads'].add(gaf.query_name)

                alignment_length = gaf.query_end - gaf.query_start
                stats['mapped_bases'] += alignment_length
                stats['total_query_bases'] += gaf.query_length
                stats['alignment_lengths'].append(alignment_length)
                stats['query_coverage'].append(alignment_length / gaf.query_length)

                if gaf.is_primary:
                    stats['primary_alignments'] += 1
                else:
                    stats['secondary_alignments'] += 1

                stats['mapq_sum'] += gaf.mapping_quality

                if gaf.strand == '+':
                    stats['strand_counts']['forward'] += 1
                else:
                    stats['strand_counts']['reverse'] += 1

                if len(gaf.path) > 1:
                    stats['multi_path_alignments'] += 1

                read_alignments[gaf.query_name].append(gaf)

    except Exception as e:
        print(f"Error processing GAF file: {e}", file=sys.stderr)
        return None

    for read, alignments in read_alignments.items():
        if len(alignments) > 1:
            stats['reads_with_multiple_alignments'] += 1

    stats['avg_mapping_quality'] = stats['mapq_sum'] / stats['total_alignments'] if stats['total_alignments'] > 0 else 0
    stats['avg_alignment_length'] = np.mean(stats['alignment_lengths']) if stats['alignment_lengths'] else 0
    stats['median_alignment_length'] = np.median(stats['alignment_lengths']) if stats['alignment_lengths'] else 0
    stats['avg_query_coverage'] = np.mean(stats['query_coverage']) if stats['query_coverage'] else 0
    stats['total_unique_reads'] = len(stats['unique_reads'])
    stats['unique_reads'] = list(stats['unique_reads'])[:5]  # Keep a few examples only

    return stats


def print_gaf_stats(stats):
    """Print formatted GAF alignment statistics"""
    print("\n===== GAF Alignment Statistics =====")
    print(f"Total alignments: {stats['total_alignments']}")
    print(f"Primary alignments: {stats['primary_alignments']} ({stats['primary_alignments'] / stats['total_alignments'] * 100:.2f}%)")
    print(f"Secondary alignments: {stats['secondary_alignments']} ({stats['secondary_alignments'] / stats['total_alignments'] * 100:.2f}%)")
    print(f"Number of uniquely aligned reads: {stats['total_unique_reads']}")
    print(f"Reads with multiple alignments: {stats['reads_with_multiple_alignments']} ({stats['reads_with_multiple_alignments'] / stats['total_unique_reads'] * 100:.2f}%)")
    print(f"Total mapped bases: {stats['mapped_bases']}")
    print(f"Total query bases: {stats['total_query_bases']}")
    print(f"Query alignment rate: {stats['mapped_bases'] / stats['total_query_bases'] * 100:.2f}%")
    print(f"Average mapping quality: {stats['avg_mapping_quality']:.2f}")
    print(f"Average alignment length: {stats['avg_alignment_length']:.2f}")
    print(f"Median alignment length: {stats['median_alignment_length']:.2f}")
    print(f"Average query coverage: {stats['avg_query_coverage'] * 100:.2f}%")
    print(f"Alignments spanning multiple paths: {stats['multi_path_alignments']} ({stats['multi_path_alignments'] / stats['total_alignments'] * 100:.2f}%)")


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_gaf.py <gaf_file> [gfa_file]")
        return

    gaf_file = sys.argv[1]
    gfa_node_dict = {}

    if len(sys.argv) > 2:
        gfa_file = sys.argv[2]
        # gfa_node_dict = parse_gfa_file(gfa_file)

    stats = analyze_gaf_file(gaf_file, gfa_node_dict)
    if stats:
        print_gaf_stats(stats)
    else:
        print("Failed to analyze GAF file.")


if __name__ == "__main__":
    main()
