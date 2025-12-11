import pyabpoa
from collections import defaultdict
import mappy as mp
from sklearn.cluster import KMeans
import warnings
import numpy as np
from sklearn.exceptions import ConvergenceWarning

from svpg.output_vcf import Candidate

def _extract_fragment_from_signature(sig, span_start, span_end, svtype):
    seq = sig.read_seq
    offset_start = span_start - sig.start
    offset_end = span_end - sig.start if svtype == 'INS' else span_end - sig.end
    left = max(0, sig.pos_read + offset_start)
    right = min(len(seq), sig.pos_read + sig.svlen + offset_end) if svtype == 'INS' else min(len(seq), sig.pos_read + offset_end)
    # left = max(0, sig.pos_read + offset_start-500)
    # right = min(len(seq), sig.pos_read + sig.svlen + offset_end+500) if svtype == 'INS' else min(len(seq), sig.pos_read + offset_end)
    middle_seq = seq[left:right]

    fragment = middle_seq
    return fragment

def _msa_consensus_for_cluster(cluster_seqs, max_cons=1, aligner=None):
    if not cluster_seqs:
        return None
    if len(cluster_seqs) == 1:
        return cluster_seqs

    if aligner is None:
        aligner = pyabpoa.msa_aligner()

    parts = sorted([(len(s), s, f'seq{i}') for i, s in enumerate(cluster_seqs)], reverse=True)
    _, seqs, names = zip(*parts)

    aln_result = aligner.msa(list(seqs), out_msa=True, out_cons=True, max_n_cons=max_cons)
    return aln_result.cons_seq

def align_with_mappy(ref_seq, query_seq, read_type):
    preset = "map-hifi" if read_type == "hifi" else "map-ont"
    aligner = mp.Aligner(seq=ref_seq,  preset=preset, fn_idx_in=None)
    if not aligner:
        raise Exception(f"Failed to load index {ref_seq}")
    aln1, aln2 = [], []
    for seq1, seq2 in query_seq:
        if not seq1:
            aln1.append(None)
        else:
            aln_primary = None
            for hit in aligner.map(seq1):
                if hit.is_primary:
                    aln_primary = (seq1, hit)
                    break
            aln1.append(aln_primary)

        if not seq2:
            aln2.append(None)
        else:
            aln_primary = None
            for hit in aligner.map(seq2):
                if hit.is_primary:
                    aln_primary = (seq2, hit)
                    break
            aln2.append(aln_primary)

    return aln1, aln2

def candidate_sv(sv, contig, ref_seq, noseqs, genotype="1/1"):
    start = sv["pos"]
    svlen = abs(sv["svlen"])
    end = start + svlen
    svtype = sv["svtype"]
    # if min_sv_length <= svlen <= max_sv_length:
    if not noseqs:
        try:
            ref_base = ref_seq[max(start - 1, 0)]
            ref_seq = ref_base if svtype == "INS" else ref_seq[max(start - 1, 0):end]
        except IndexError:
            return
        alt_seq = ref_base if svtype == "DEL" else ref_base + sv['seq']
    else:
        ref_seq = "N"
        alt_seq = "<INS>" if svtype == "INS" else "<DEL>"

    return Candidate(contig, start, end, svtype, [], ref_seq=ref_seq, alt_seq=alt_seq, genotype=genotype)

def phase_reads_by_cluster_similarity(current_signature_clusters):
    """
    phase reads into 2 haps based on their signature cluster similarity using k-means clustering.
    """
    # Step 1: construct read_name -> sigs mapping
    read_to_sigs = defaultdict(list)
    for cluster in current_signature_clusters:
        for sig in cluster:
            read_to_sigs[sig.read_name].append(sig)

    read_names = list(read_to_sigs.keys())
    n_reads = len(read_names)
    n_clusters = len(current_signature_clusters)

    # Step 2: construct feature matrix （n_reads x n_clusters）
    # sv of read i in cluster j -> feature_matrix[i][j] = 1 else 0
    feature_matrix = np.zeros((n_reads, n_clusters), dtype=int)

    cluster_read_names = [set(sig.read_name for sig in cluster)
                          for cluster in current_signature_clusters]

    for i, read_name in enumerate(read_names):
        for j, cluster_reads in enumerate(cluster_read_names):
            if read_name in cluster_reads:
                feature_matrix[i, j] = 1

    # Step 3: k-means
    if n_reads >= 2:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            kmeans = KMeans(n_clusters=2, random_state=42, n_init=20)
            labels = kmeans.fit_predict(feature_matrix)
    else:
        labels = np.zeros(n_reads, dtype=int)

    # Step 4: read_name -> hap mapping
    read_to_hap = {read: (label + 1) for read, label in zip(read_names, labels)}  # hap1 / hap2

    # Step 5: construct hap -> seqs mapping
    hap_to_seqs = defaultdict(list)
    for read_name, hap_id in read_to_hap.items():
        sigs = read_to_sigs[read_name]
        span_start = max(min(s.pos_read for s in sigs)-2000, 0)
        span_end = min(max(s.pos_read for s in sigs)+2000, len(sigs[-1].read_seq))
        if sigs[-1].type == 'INS':
            frag_seq = sigs[-1].read_seq[span_start: span_end+sigs[-1].svlen]
        else:
            frag_seq = sigs[-1].read_seq[span_start: span_end]
        hap_to_seqs[hap_id].append(frag_seq)

    return list(hap_to_seqs.values())


def _extract_sv_from_alignment(aln):
    a1_seq, a1_res = aln[0], aln[1]
    if a1_res.mapq < 20:
        return []
    res = []
    ref_pos = a1_res.r_st  # 0-based
    read_pos = 0
    # CIGAR tuples: (operation, length)
    cigartuples = [(length, op) for length, op in a1_res.cigar]
    for length, op in cigartuples:
        # 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
        if op in (0, 7, 8):
            ref_pos += length
            read_pos += length
        elif op == 1 and length >= 50:  # INS in read (relative to reference)
            # extract inserted sequence from read
            insersion_seq = a1_seq[read_pos:read_pos + length]
            res.append({
                        'pos': ref_pos,  # insertion is between ref_pos-1 and ref_pos
                        'svtype': 'INS',
                        'svlen': length,
                        'read_pos': read_pos,
                        'seq': insersion_seq
                        })
            read_pos += length
        elif op == 2 and length >= 50:  # DEL from reference
            res.append({
                        'pos': ref_pos,
                        'svtype': 'DEL',
                        'svlen': length,
                        'read_pos': read_pos,
                        })
            ref_pos += length
        elif op in (4, 5):  # soft/hard clip
            read_pos += length
        else:
            # skip others
            pass
    return res

def run_align(merged_intervals, ref_seq, options):
    haps, beds, beds_sv = [], [], []
    # complex regions realignment
    for contig, start, end, svtype, current_signature_clusters, _ in merged_intervals:
        if len(current_signature_clusters) < 3:
            continue
        hap_to_seqs = phase_reads_by_cluster_similarity(current_signature_clusters)

        if not hap_to_seqs:
            continue
        cons_hp1 = _msa_consensus_for_cluster(hap_to_seqs[0])[0]
        if len(hap_to_seqs) == 1:
            continue
        cons_hp2 = _msa_consensus_for_cluster(hap_to_seqs[1])[0]
        haps.append((cons_hp1, cons_hp2))
        beds.append((start, end))
        beds_sv.append(current_signature_clusters)

    sv_candidates = []
    min_sv_length, noseqs = options.min_sv_size, options.noseq
    aln1, aln2 = align_with_mappy(ref_seq, haps, options.read)
    for a1, a2, bed, align_sv in zip(aln1, aln2, beds, beds_sv):
        tmp_sv = []
        start, end = bed

        svs1 = _extract_sv_from_alignment(a1) if a1 else []
        svs2 = _extract_sv_from_alignment(a2) if a2 else []

        # ---- filter svs within region ±1kb ----
        svs1 = [sv for sv in svs1 if start-1000 <= sv['pos'] <= end+1000]
        svs2 = [sv for sv in svs2 if start-1000 <= sv['pos'] <= end+1000]

        # Case A: both haplotypes have SVs
        if svs1 and svs2:
            for sv1 in svs1:
                matched = False
                for sv2 in svs2:
                    if sv1['svtype'] == sv2['svtype']:
                        len1, len2 = abs(sv1['svlen']), abs(sv2['svlen'])
                        if min(len1, len2) / max(len1, len2) > 0.7:
                            # hp1 match hp2 → 1/1
                            sv_candidates.append(candidate_sv(
                                sv1, contig, ref_seq, noseqs, genotype="1/1"
                            ))
                            matched = True
                            svs2.remove(sv2)
                            break

                if not matched:
                    # hp1 not match hp2 → 0/1
                    sv_candidates.append(candidate_sv(
                        sv1, contig, ref_seq, noseqs, genotype="0/1"
                    ))

            # hp2 not match hp1 → 0/1
            for sv2 in svs2:
                sv_candidates.append(candidate_sv(
                    sv2, contig, ref_seq, noseqs, genotype="0/1"
                ))

        # Case B: only hp1 has SVs → 0/1
        elif svs1:
            for sv1 in svs1:
                sv_candidates.append(candidate_sv(
                    sv1, contig, ref_seq, noseqs, genotype="0/1"
                ))

        # Case C:  only hp2 has SVs → 0/1
        elif svs2:
            for sv2 in svs2:
                sv_candidates.append(candidate_sv(
                    sv2, contig, ref_seq, noseqs, genotype="0/1"
                ))

    return sv_candidates