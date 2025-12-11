from random import sample

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

def form_bins(sv_signatures, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_signatures = sorted(sv_signatures, key=lambda evi: evi.get_key())
    grouped_bin, current_group = [], []
    for i, sv_sig in enumerate(sorted_signatures):
        if i == 0:
            current_group.append(sv_sig)
        else:
            if current_group[-1].downstream_distance_to(sv_sig) < max_distance:
                current_group.append(sv_sig)
            else:
                grouped_bin.append(current_group)
                current_group = [sv_sig]
    if current_group:
        grouped_bin.append(current_group)

    bin_depth = [len(bin) for bin in grouped_bin]
    try:
        mean_depth = sum(bin_depth) / len(bin_depth)
    except ZeroDivisionError:
        mean_depth = 0

    return grouped_bin, mean_depth

def span_position_distance(signature1, signature2, type, param):
    span1 = signature1.get_source()[2] - signature1.get_source()[1]
    span2 = signature2.get_source()[2] - signature2.get_source()[1]
    center1 = (signature1.get_source()[1] + signature1.get_source()[2]) // 2
    center2 = (signature2.get_source()[1] + signature2.get_source()[2]) // 2
    position_distance = abs(center1 - center2) // (max(span1, span2) * param)
    node_distance = 0
    if type == 'BND':
        if signature1.get_destination()[0] != signature2.get_destination()[0]:
            return 999
        dir1 = (signature1.source_direction, signature1.dest_direction)
        dir2 = (signature2.source_direction, signature2.dest_direction)
        if dir1 != dir2:
            return 999
        dist1 = abs(signature1.get_source()[1] - signature2.get_source()[1])
        dist2 = abs(signature1.get_destination()[1] - signature2.get_destination()[1])
        if max(dist1, dist2) == 0:
            return 0

        return dist1 // (max(dist1, dist2) * param) + dist2 // (max(dist1, dist2) * param)
    elif type == 'INS' or type == 'DEL':
        if signature1.node_ls and signature2.node_ls:
            path1, path2 = signature1.node_ls, signature2.node_ls
            node_distance = 0.05 * (1 - len(set(path1).intersection(set(path2))) / len(set(path1).union(set(path2))))
            # node_distance = len(set(path1).symmetric_difference(set(path2))) / max(len(path1), len(path2)) * 0.1

    span_distance = abs(span1 - span2) / max(span1, span2)

    return position_distance + span_distance + node_distance

def cluster_data(bins, mean_len):
    clusters_final = []
    for bin in bins:
        if len(bin) == 1:
            clusters_final.append(bin)
            continue
        elif len(bin) > 100:
            partition_sample = sample(bin, 100)
        else:
            partition_sample = bin

        element_type = partition_sample[0].type
        param = (abs(len(partition_sample) - mean_len)) / max(len(partition_sample), mean_len)+mean_len

        distance_data = []
        for i in range(len(partition_sample) - 1):
            for j in range(i + 1, len(partition_sample)):
                distance_data.append(
                    span_position_distance(partition_sample[i], partition_sample[j], element_type, param))
        Z = linkage(np.array(distance_data), method="average")
        cluster_indices = list(fcluster(Z, 0.3, criterion='distance'))
        new_clusters = [[] for i in range(max(cluster_indices))]
        for signature_index, cluster_index in enumerate(cluster_indices):
            new_clusters[cluster_index - 1].append(partition_sample[signature_index])
        clusters_final.extend(new_clusters)

    return clusters_final
