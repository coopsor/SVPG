from math import log10
import numpy as np
# import pysam

err = 0.2
prior = float(1/3)
Genotype = ["0/0", "0/1", "1/1"]

def threshold_ref_count(num):
    if num <= 2:
        return 20*num
    elif 3 <= num <= 5:
        return 9*num
    elif 6 <= num <= 15:
        return 7*num
    else:
        return 5*num

def log10sumexp(log10_probs):
    # Normalization of Genotype likelihoods
    m = max(log10_probs)
    return m + log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
    # Adjust the Genotype likelihoods
    log10_probs = np.array(log10_probs)
    lse = log10sumexp(log10_probs)
    return np.minimum(log10_probs - lse, 0.0)

def rescale_read_counts(c0, c1, max_allowed_reads=100):
    """Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
    Total = c0 + c1
    if Total > max_allowed_reads:
        c0 = int(max_allowed_reads * float(c0/Total))
        c1 = max_allowed_reads - c0
    return c0, c1

def cal_GL(c0, c1, type, platform):
    if platform == "hifi":
        err = 0.05 if type == "INS" else 0.05
    else:
        err = 0.1 if type == "INS" else 0.1
    # Approximate adjustment of events with larger read depth
    c0, c1 = rescale_read_counts(c0, c1)

    ori_GL00 = np.float64(pow((1 - err), c0) * pow(err, c1) * (1 - prior) / 2)
    ori_GL11 = np.float64(pow(err, c0) * pow((1 - err), c1) * (1 - prior) / 2)
    ori_GL01 = np.float64(pow(0.5, c0 + c1) * prior)

    # normalized genotype likelihood
    prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10 * log10(i))) for i in GL_P]
    GQ = [int(-10 * log10(GL_P[1] + GL_P[2])), int(-10 * log10(GL_P[0] + GL_P[2])), int(-10 * log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10 * log10(GL_P[0]), 1))

    return Genotype[prob.index(max(prob))], "%d,%d,%d" % (PL[0], PL[1], PL[2]), max(GQ), QUAL

def genotype(candidates, type, options):
    # bam = pysam.AlignmentFile(options.bam, threads=options.num_threads)
    from svpg.main import bam
    for nr, candidate in enumerate(candidates):
        reads_supporting_variant = set(candidate.members)
        max_bias = 1000
        contig, start, end = candidate.get_source()
        contig_length = bam.get_reference_length(contig)
        try:
            alignment_it = bam.fetch(contig=contig, start=max(0, start-max_bias), stop=min(contig_length, end+max_bias))
        except ValueError:
            continue
        #  Count reads that overlap the locus and therefore support the reference
        aln_no = 0
        reads_supporting_reference = set()

        while aln_no < 500:
            try:
                current_alignment = next(alignment_it)
            except StopIteration:
                break
            if current_alignment.query_name in reads_supporting_variant:
                continue
            if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
                continue
            aln_no += 1
            if type == "DEL":
                minimum_overlap = min((end - start) / 2, 2000)
                if (current_alignment.reference_start < (end - minimum_overlap) and current_alignment.reference_end > (end + max_bias) or
                    current_alignment.reference_start < (start - max_bias) and current_alignment.reference_end > (start + minimum_overlap)):
                    reads_supporting_reference.add(current_alignment.query_name)
            else:
                if current_alignment.reference_start < (start - max_bias) and current_alignment.reference_end > (end + max_bias):
                    reads_supporting_reference.add(current_alignment.query_name)

        GT, GL, GQ, QUAL = cal_GL(len(reads_supporting_reference), len(reads_supporting_variant), type, options.read)

        candidate.support_fraction = len(reads_supporting_variant) / (
                len(reads_supporting_variant) + len(reads_supporting_reference))
        candidate.genotype = GT
        candidate.ref_reads = len(reads_supporting_reference)
        candidate.alt_reads = len(reads_supporting_variant)

    return candidates