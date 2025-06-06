import os
from collections import defaultdict

ref_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y',
            'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
            'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

linear_regions = defaultdict(list)
gfa_edges = defaultdict(dict)

# Parse GFA file to extract linear regions and edges
with open('/data/huheng/dataset/human/pangenome/GRCh37-91c.r559.gfa', "r") as fp:
    for line in fp:
        tokens = line.strip().split('\t')
        if tokens[0] == 'S':  # Segment line
            if int(tokens[6][5:]) == 0:  # Only include reference path
                linear_regions[tokens[4][5:]].append([int(tokens[5][5:]), int(tokens[5][5:]) + int(tokens[3][5:])])
        if tokens[0] == 'L':  # Link line
            gfa_edges.setdefault(tokens[1], []).append(tokens[3])

# Compute "depth" of nodes in the region (i.e., branching complexity)
def calculate_depth(gfa_edges, start_node, end_node):
    node_depth = 0
    for index in range(start_node, end_node):
        node = 's' + str(index + 1)
        if node not in gfa_edges:
            continue
        node_depth += len(gfa_edges[node])
    return node_depth

# Binary search to find the first overlapping subregion
def find_start_subregion(subregions, query_start):
    low, high = 0, len(subregions) - 1
    while low <= high:
        mid = (low + high) // 2
        if subregions[mid][1] < query_start:
            low = mid + 1
        else:
            high = mid - 1
    return low

# Merge overlapping subregions with the query region
def merge_subregions(subregions, query_start, query_end):
    i = find_start_subregion(subregions, query_start)
    if i == len(subregions):
        return None
    merged_start, merged_end = subregions[i][0], subregions[i][1]
    j = i
    while j < len(subregions) and subregions[j][0] <= query_end:
        merged_end = max(merged_end, subregions[j][1])
        j += 1
    if query_end > merged_end:
        merged_end = query_end
    return (merged_start, merged_end, i, j)

# Output depth of each 5kb region and mark low-depth regions separately
with open('node_depth.bed', 'w') as f:
    with open('node_low_depth.bed', 'w') as f1:
        node_depth_list = []
        for chrom, chrom_list in linear_regions.items():
            if chrom not in ref_list:
                continue
            chrom_length = chrom_list[-1][1]
            for start_coord in range(0, chrom_length, 5000):
                end_coord = start_coord + 5000
                if end_coord > chrom_length:
                    end_coord = chrom_length
                result = merge_subregions(chrom_list, start_coord, end_coord)
                if result is not None:
                    node_start_coord, node_end_coord, node_start_index, node_end_index = result
                else:
                    break
                node_depth = calculate_depth(gfa_edges, node_start_index, node_end_index)
                node_depth_list.append([start_coord, end_coord, node_depth])
            for rec in node_depth_list:
                start_coord, end_coord, node_depth = rec
                if node_depth >= 4:
                    f.write('{0}\t{1}\t{2}\t{3}\n'.format(str(chrom), start_coord, end_coord, node_depth))
                else:
                    f1.write('{0}\t{1}\t{2}\t{3}\n'.format(str(chrom), start_coord, end_coord, node_depth))

# Merge complex graph regions with depth >= 4, allowing 1000 bp gap
os.system('sort -k1,1 -k2,2n node_depth.bed | bedtools merge -d 1000 -c 4 -o sum > complex_graph.bed')
