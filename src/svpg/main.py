import pickle
import re
import sys
import os
import logging
import time
from multiprocessing import Pool
import pysam
import numpy as np
from time import strftime, localtime
import subprocess

from svpg.input_parsing import parse_arguments
from svpg.SVCollect import read_bam
from svpg.SVCluster import form_bins, cluster_data
from svpg.SVPan import read_gaf, read_gaf_pan
from svpg.util import read_gfa, find_sequence_file
from svpg.output_vcf import consolidate_clusters_unilocal, write_final_vcf
from svpg.SVGenotype import genotype
from svpg.graph_augment import augment_pipe
from svpg.realign import run_align

options = parse_arguments()
ref_genome = pysam.FastaFile(options.ref)


def multi_process(total_len, step, args=None):
    num_threads = min(options.num_threads, max(1, total_len // 100))

    chunk_size = total_len // num_threads
    chunks = []

    for i in range(num_threads):
        start = i * chunk_size
        end = start + chunk_size if i < num_threads - 1 else total_len
        if step == 'read_bam':
            chunks.append((args, start, end, options))
        elif step == 'cluster':
            chunks.append((args[0][start:end], args[1]))
        else:
            chunks.append((args[0][start:end], args[1], options))

    with Pool(processes=num_threads) as pool:
        if step == 'read_bam':
            results = pool.starmap(read_bam, chunks)
        elif step == 'read_gaf':
            results = pool.starmap(read_gaf, chunks)
        elif step == 'read_gaf_pan':
            results = pool.starmap(read_gaf_pan, chunks)
        elif step == 'realign':
            results = pool.starmap(run_align, chunks)
        elif step == 'cluster':
            results = pool.starmap(cluster_data, chunks)
        else:
            results = pool.starmap(genotype, chunks)

    return [item for sublist in results for item in sublist]


def read_in_chunks(file_object, chunk_size=102400):
    while True:
        lines = []
        for _ in range(chunk_size):
            line = file_object.readline().decode('utf-8').strip()
            if not line:
                break
            lines.append(line)
        if not lines:
            break
        yield lines

def recall_task(positions, adjacent, signature_clusters):
    merged_intervals = []  # [(chrom, start, end, svtype, [cluster_idx,...])]
    current_start = current_end = current_contig = current_svtype = None
    current_indices = []

    ends = [np.median([m.end for m in c]) for c in signature_clusters]
    for i in range(len(positions)):
        if not adjacent[i]:
            continue

        pos = positions[i]
        contig = signature_clusters[i][0].contig
        svtype = signature_clusters[i][0].type
        end_pos = ends[i] if svtype != "INS" else pos

        if current_start is None:
            current_contig = contig
            current_svtype = svtype
            current_start = pos
            current_end = end_pos
            current_indices = [i]
            current_signature_clusters = [signature_clusters[i]]
            continue

        last_i = current_indices[-1]
        last_pos = positions[last_i]

        if contig == current_contig and abs(pos - last_pos) < 1000:
            # expand current interval
            current_indices.append(i)
            current_signature_clusters.append(signature_clusters[i])
            current_end = end_pos  # update end position
        else:
            # save current interval and start a new one
            merged_intervals.append((current_contig, int(current_start), int(current_end), current_svtype, current_signature_clusters, current_indices))

            current_contig = contig
            current_svtype = svtype
            current_start = pos
            current_end = end_pos
            current_indices = [i]
            current_signature_clusters = [signature_clusters[i]]

    if current_start is not None:
        merged_intervals.append(
            (current_contig, int(current_start), int(current_end),
             current_svtype, current_signature_clusters, current_indices)
        )

    chrom_merged = {}
    uncalled_clusters, recalled_sv = [], []
    for contig, start, end, svtype, sigs, idx_list in merged_intervals:
        chrom_merged.setdefault(contig, []).append((contig, start, end, svtype, sigs, idx_list))

    for chrom, intervals in chrom_merged.items():
        ref_seq = ref_genome.fetch(chrom)
        recall_candidates = multi_process(len(intervals), 'realign', (intervals, ref_seq))
        seen = set()

        j = 0
        n = len(recall_candidates)

        for contig, start, end, svtype, sigs, idx_list in intervals:
            found = False
            contained_svs = []

            while j < n and recall_candidates[j].end < start:
                j += 1

            k = j
            while k < n and recall_candidates[k].start <= end + 1000:
                sv = recall_candidates[k]
                if start - 1000 <= sv.start <= end + 1000:
                    found = True
                    contained_svs.append(sv)
                k += 1

            if not found:
                for i in idx_list:
                    uncalled_clusters.append(signature_clusters[i])
            else:
                for sv in contained_svs:
                    sv_id = (sv.start, sv.end)
                    if sv_id not in seen:
                        recalled_sv.append(sv)
                        seen.add(sv_id)

    logging.info(f"Recalled {len(recalled_sv)} SVs and {len(uncalled_clusters)} uncalled clusters.")
    return recalled_sv, uncalled_clusters

def main():
    # Set up logging
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Ensure the base directory exists
    os.makedirs(options.working_dir, exist_ok=True)

    # Create log file
    fileHandler = logging.FileHandler(
        "{0}/SVPG_{1}.log".format(options.working_dir, strftime("%y%m%d_%H%M%S", localtime())), mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    logging.info("******************************** Start SVPG ******************************")

    logging.info("CMD: python3 {0}".format(" ".join(sys.argv)))
    logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))

    gfa_node = read_gfa(options.gfa)

    if options.contigs is None:
        options.contigs = [ctg for ctg in ref_genome.references if re.match(r'^(chr)?[0-9XYM]+$', ctg)]

    if options.max_merge_threshold is None:
        if options.read == 'hifi':
            options.max_merge_threshold = 50
        else:
            options.max_merge_threshold = 500

    for arg in vars(options):
        logging.info("PARAMETER: {0}, VALUE: {1}".format(arg, getattr(options, arg)))

    pan_signatures = []
    if options.sub == 'call':
        logging.info("MODE: call")
        logging.info("INPUT: {0}".format(os.path.abspath(options.bam)))
        logging.info("***************** Collect SV signatures *****************")

        try:
            bam = pysam.AlignmentFile(options.bam, threads=options.num_threads)
            bam.check_index()
        except ValueError:
            logging.warning(
                "Input BAM file is missing a valid index. Please generate with 'samtools faidx'.")
        except AttributeError:
            logging.warning(
                "pysam's check_index raised an Attribute error. Something is wrong with the input BAM file.")
            return

        bam_signatures = []
        ref_list = bam.get_index_statistics()
        for ref in ref_list:
            if ref.mapped == 0:
                continue
            if ref[0] in options.contigs:
                ref_len = bam.get_reference_length(ref[0])
                logging.info("Processing ref {0}...".format(ref[0]))
                bam_signatures.extend(multi_process(ref_len, 'read_bam', ref[0]))
                logging.info("Processed ref {0}...".format(ref[0]))

        logging.info("****************************** Graph Mapping ******************************")

        # with open(options.working_dir + '/sv_signatures.pkl', 'wb') as temp:
        #     pickle.dump(bam_signatures, temp)
        # with open(options.working_dir + '/sv_signatures.pkl', 'rb') as f:
        #     bam_signatures = pickle.load(f)

        deletion_signatures = [ev for ev in bam_signatures if ev.type == "DEL"]
        insertion_signatures = [ev for ev in bam_signatures if ev.type == "INS"]
        # logging.info("Found {0} signatures for deleted regions.".format(len(deletion_signatures)))
        # logging.info("Found {0} signatures for inserted regions.".format(len(insertion_signatures)))

        signature_clusters = []
        for element_signature in [insertion_signatures, deletion_signatures]:
            if not element_signature:
                continue
            signature_bin, bin_depth = form_bins(element_signature, 1000)
            if bin_depth == 0:
                logging.warning("No signatures found in the current bin. Skipping clustering for this bin.")
                continue
            signature_clusters.extend(multi_process(len(signature_bin), 'cluster', (signature_bin, bin_depth)))

        signature_clusters = sorted(signature_clusters, key=lambda x: (x[0].contig, np.median([m.start for m in x])))
        positions = [np.median([m.start for m in c]) for c in signature_clusters]

        n = len(positions)
        adjacent = [False] * n

        for i in range(n):
            if i > 0 and abs(positions[i] - positions[i - 1]) < 1000:
                adjacent[i] = adjacent[i - 1] = True
            if i < n - 1 and abs(positions[i + 1] - positions[i]) < 1000:
                adjacent[i] = adjacent[i + 1] = True
        close_clusters = [signature_clusters[i] for i in range(n) if adjacent[i]]

        if options.realign:
            logging.info("Realignment enabled: Merging adjacent clusters for realignment.")
            recalled_sv, uncalled_clusters = recall_task(positions, adjacent, signature_clusters)

        refine_bins = [signature_clusters[i] for i in range(n) if not adjacent[i]]
        refine_sigs = [sig for group in refine_bins for sig in group]
        sig_read = 'signatures'

        fasta_file = open(options.working_dir + f'/{sig_read}.fa', 'w')
        for sig_index, sig in enumerate(refine_sigs):
            # adjac_distance = max(min(5000, sig.svlen*3), 2000)
            adjac_distance = 2000
            if sig.signature == 'suppl':
                read_seq = sig.read_seq
            else:
                if sig.pos_read < adjac_distance:
                    read_seq = sig.read_seq[0:sig.pos_read + sig.svlen + adjac_distance]
                else:
                    read_seq = sig.read_seq[
                               sig.pos_read - adjac_distance:sig.pos_read + sig.svlen + adjac_distance] if sig.pos_read else sig.read_seq

            ref_suppl1, ref_suppl2 = '', ''
            svtype = sig.type
            if sig.pos_read < adjac_distance:
                try:
                    ref_suppl1 = ref_genome.fetch(sig.contig, sig.start - adjac_distance, sig.start - sig.pos_read)
                except ValueError:
                    ref_suppl1 = ''
            if svtype == 'DEL' and sig.pos_read + adjac_distance > len(sig.read_seq):
                try:
                    ref_suppl2 = ref_genome.fetch(sig.contig, sig.end + (len(sig.read_seq) - sig.pos_read),
                                                  sig.end + adjac_distance)
                except ValueError:
                    ref_suppl2 = ref_genome.fetch(sig.contig, sig.end + (len(sig.read_seq) - sig.pos_read),
                                                  len(ref_genome.fetch(sig.contig)))
            elif svtype == 'INS' and sig.pos_read + sig.svlen + adjac_distance > len(sig.read_seq):
                try:
                    ref_suppl2 = ref_genome.fetch(sig.contig,
                                                  sig.start + len(sig.read_seq) - sig.pos_read - sig.svlen,
                                                  sig.start + adjac_distance)
                except ValueError:
                    ref_suppl2 = ref_genome.fetch(sig.contig,
                                                  sig.start + len(sig.read_seq) - sig.pos_read - sig.svlen,
                                                  len(ref_genome.fetch(sig.contig)))

            read_seq = ref_suppl1 + read_seq + ref_suppl2
            pos_ref = str(sig.contig) + ':' + str(sig.start) + ':' + str(sig.end)
            read_info = f"{sig.read_name}@{svtype}@{pos_ref}@{sig.alt_seq}" if svtype == 'INS' else f"{sig.read_name}@{svtype}@{pos_ref}"
            fasta_file.write(f'>{read_info}\n{read_seq}\n')

        fasta_file.close()

        if options.read == 'hifi':
            os.system(
                f'minigraph -t {options.num_threads} -cx asm --vc --secondary yes {options.gfa} {options.working_dir}/{sig_read}.fa > {options.working_dir}/{sig_read}.gaf')
        else:
            os.system(
                f'minigraph -t {options.num_threads} -cx lr --vc --secondary yes {options.gfa} {options.working_dir}/{sig_read}.fa > {options.working_dir}/{sig_read}.gaf')

        logging.info("*************** Collect signatures from pangenome-reference ***************")

        pan_signatures = read_gaf(gfa_node, options)
        # with open(options.working_dir + f'/{sig_read}.gaf', 'rb') as f:
        #     for chunk_index, lines in enumerate(read_in_chunks(f, chunk_size=200000000)):
        #         logging.info(f"Processing chunks {chunk_index + 1}")
        #         pan_signatures.extend(multi_process(len(lines), 'read_gaf', (lines, gfa_node)))
        #         logging.info(f"Processed chunks {chunk_index + 1}")

    elif options.sub == 'graph-call':
        logging.info("MODE: graph-call")
        logging.info("INPUT: {0}".format(os.path.abspath(options.gaf)))
        logging.info("*************** Collect SV signatures from pangenome ***************")

        pan_signatures = read_gaf_pan(gfa_node, options)
        # with open(options.gaf, 'rb') as f:
        #     for chunk_index, lines in enumerate(read_in_chunks(f, chunk_size=200000000)):
        #         logging.info(f"Processing chunk {chunk_index + 1}")
        #         pan_signatures.extend(multi_process(len(lines), 'read_gaf_pan', (lines, gfa_node)))
        #         logging.info(f"Processed chunks {chunk_index + 1}")
        sig_read = 'signatures_test'
    elif options.sub == 'augment':
        logging.info("MODE: augment")
        logging.info("*************** Collect SVs from pangenome ***************")
        start_time = time.time()

        base_dir = options.working_dir
        if not options.skip_call:
            filelist_path = os.path.join(base_dir, "filelist.tsv")
            if os.path.exists(filelist_path):
                os.remove(filelist_path)

            sample_paths_to_process = []
            if options.sample_list:
                try:
                    with open(options.sample_list, 'r') as f:
                        sample_paths_to_process = [line.strip() for line in f if line.strip()]
                except FileNotFoundError:
                    logging.error(f"Sample list file not found: {options.sample_list}")
                    raise RuntimeError(f"Sample list file not found: {options.sample_list}")
            else:
                # Fallback to directory scanning if sample_list is not provided
                for entry in os.scandir(base_dir):
                    if entry.is_dir() and entry.name.startswith("sample"):
                        # Assuming FASTA file is directly inside the sample directory and named as {prefix}.fasta
                        file_type = find_sequence_file(entry)
                        if not file_type:
                            logging.warning(f"Expected FASTA file {file_type} not found. Skipping {entry.name}.")
                            raise RuntimeError(f"FASTA file not found for sample: {entry.name}")

                        fasta_path_in_dir = os.path.join(entry.path, f"{entry.name}{file_type}")
                        sample_paths_to_process.append(fasta_path_in_dir)
            if not sample_paths_to_process:
                logging.error("No sample paths to process. Please check your sample list or directory structure.")
                raise RuntimeError("No sample paths to process.")
            else:
                logging.info(f"Found {len(sample_paths_to_process)} samples to process.")
            with open(filelist_path, "a") as filelist:
                for fasta_file_path in sample_paths_to_process:
                    try:
                        # Get the directory of the fasta file and its prefix
                        sample_dir = os.path.dirname(fasta_file_path)
                        prefix = os.path.basename(sample_dir) if sample_dir else \
                            os.path.splitext(os.path.basename(fasta_file_path))[0]

                        original_cwd = os.getcwd()
                        os.chdir(sample_dir)

                        fasta_file_name = os.path.basename(fasta_file_path)

                        logging.info(f"Start call SVs from {prefix}")
                        file_size = os.path.getsize(fasta_file_name)
                        coverage = file_size // (1024 * 1024 * 1024) // 3.1
                        hifi_support_map = [
                            (0, 5, 1), (5, 15, 2), (15, 25, 3), (25, 50, 4), (50, float("inf"), 5)
                        ]
                        ont_support_map = [
                            (0, 5, 2), (5, 15, 3), (15, 25, 4), (25, 50, 5), (50, float("inf"), 10)
                        ]
                        support_map = hifi_support_map if options.read == 'hifi' else ont_support_map
                        support = next(val for low, high, val in support_map if low <= coverage <= high)

                        gaf_file = f"{prefix}.gaf"
                        if not os.path.exists(gaf_file):
                            if options.read == 'hifi':
                                cmd_align = f"minigraph -t{options.num_threads} -cxasm --vc --secondary yes {options.gfa} {fasta_file_name} > {gaf_file}"
                            else:
                                cmd_align = f"minigraph -t{options.num_threads} -cxlr --vc --secondary yes {options.gfa} {fasta_file_name} > {gaf_file}"
                            os.system(cmd_align)

                        var_file = options.vcf_out
                        cmd_call = [
                            "python", __file__, "graph-call",
                            "--read", options.read,
                            "-s", str(support),
                            "-t",str(options.num_threads),
                            "--working_dir", './',  # This should refer to the current sample directory
                            "--ref", options.ref,
                            "--gfa", options.gfa,
                            "--gaf", gaf_file,
                            "-o", var_file,
                            "--min_sv_size", str(options.min_sv_size),
                            "--max_sv_size", str(options.max_sv_size),
                            "--types", 'DEL,INS'
                        ]
                        try:
                            subprocess.run(cmd_call, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                           text=True)
                            pysam.tabix_compress(var_file, f"{var_file}.gz", force=True)
                            pysam.tabix_index(f"{var_file}.gz", preset="vcf", force=True)
                        except subprocess.CalledProcessError:
                            logging.error(
                                f"'{prefix}' encountered an error while running the SVs call.\
                                Please check the logs in the sample directory: {sample_dir}"
                            )
                            raise RuntimeError(f"Error occurred for sample: {prefix}") from e

                        # write the VCF path to filelist.tsv
                        vcf_path = os.path.join(sample_dir, f"{var_file}.gz")
                        if os.path.exists(vcf_path):
                            filelist.write(f"{vcf_path}\n")

                    except Exception as e:
                        logging.error(f"Failed to process sample from path {fasta_file_path}: {e}")
                    finally:
                        # Always change back to the original working directory
                        os.chdir(original_cwd)

        call_time = time.time()
        logging.info(f"SVs call time: {call_time - start_time:.2f} seconds")

        logging.info("*************** Augment pangenome graph ***************")
        augment_pipe(base_dir, options.ref, options.gfa, options.out)
        end_time = time.time()
        logging.info(f"Graph augment time: {end_time - call_time:.2f} seconds")
        logging.info(f"Total time: {end_time - start_time:.2f} seconds")

        return

    # with open(options.working_dir + f'/{sig_read}_test.pkl', 'wb') as temp:
    #     pickle.dump(pan_signatures, temp)
    # with open(options.working_dir + f'/{sig_read}_test.pkl', 'rb') as f:
    #     pan_signatures = pickle.load(f)

    deletion_signatures = [ev for ev in pan_signatures if ev.type == "DEL"]
    insertion_signatures = [ev for ev in pan_signatures if ev.type == "INS"]
    duplication_signatures = [ev for ev in pan_signatures if ev.type == "DUP"]
    inversion_signatures = [ev for ev in pan_signatures if ev.type == "INV"]
    breakend_signatures = [ev for ev in pan_signatures if ev.type == "BND"]

    logging.info("Found {0} signatures for deleted regions.".format(len(deletion_signatures)))
    logging.info("Found {0} signatures for inserted regions.".format(len(insertion_signatures)))
    logging.info("Found {0} signatures for duplicated regions.".format(len(duplication_signatures)))
    logging.info("Found {0} signatures for inverted regions.".format(len(inversion_signatures)))
    logging.info("Found {0} signatures for translated regions.".format(len(breakend_signatures)))

    pan_clusters = []
    for element_signature in [deletion_signatures, insertion_signatures, duplication_signatures, inversion_signatures, breakend_signatures]:
        if not element_signature:
            continue
        signature_bin, bin_depth = form_bins(element_signature, 1000)
        if bin_depth == 0:
            logging.warning("No signatures found in the current bin. Skipping clustering for this bin.")
            continue

        pan_clusters.extend(multi_process(len(signature_bin), 'cluster', (signature_bin, bin_depth)))

    if options.sub == 'call':
        if options.realign:
            pan_clusters = pan_clusters + uncalled_clusters
        else:
            pan_clusters = pan_clusters + close_clusters
    pan_clusters = [cluster for cluster in pan_clusters if len(cluster) >= options.min_support]
    chrom_results = {}
    for sig in pan_clusters:
        chrom_results.setdefault(sig[0].contig, []).append(sig)

    logging.info("********************************** SVCALL *********************************")

    sv_candidate = []
    for contig in chrom_results:
        if contig in options.contigs:
            ref_chrom_seq = ref_genome.fetch(contig)
            sv_candidate.extend(sorted(
                consolidate_clusters_unilocal(chrom_results[contig], ref_chrom_seq, options, cons=options.alt_consensus),
                key=lambda cluster: (cluster.contig, cluster.start)))

    deletion_candidates = [i for i in sv_candidate if i.type == 'DEL']
    insertion_candidates = [i for i in sv_candidate if i.type == 'INS']
    duplication_candidates = [i for i in sv_candidate if i.type == 'DUP']
    breakend_candidates = [i for i in sv_candidate if i.type == 'BND']

    logging.info("Final deletion candidates: {0}".format(len(deletion_candidates)))
    logging.info("Final insertion candidates: {0}".format(len(insertion_candidates)))
    logging.info("Final duplication candidates: {0}".format(len(duplication_candidates)))
    logging.info("Final breakend candidates: {0}".format(len(breakend_candidates)))

    if options.sub == 'call' and not options.skip_genotype:
        logging.info("********************************* GENOTYPE ********************************")
        logging.info("Genotyping deletions..")
        deletion_candidates = multi_process(len(deletion_candidates), 'genotype', (deletion_candidates, "DEL"))
        logging.info("Genotyping insertions..")
        insertion_candidates = multi_process(len(insertion_candidates), 'genotype', (insertion_candidates, "INS"))
        logging.info("Genotyping duplications..")
        duplication_candidates = multi_process(len(duplication_candidates), 'genotype', (duplication_candidates, "DUP"))
        logging.info("Genotyping breakends..")
        breakend_candidates = multi_process(len(breakend_candidates), 'genotype', (breakend_candidates, "BND"))

    if options.sub == 'call' and options.realign:
        deletion_candidates_recall = [i for i in recalled_sv if i.type == 'DEL']
        insertion_candidates_recall = [i for i in recalled_sv if i.type == 'INS']
        deletion_candidates += deletion_candidates_recall
        insertion_candidates += insertion_candidates_recall

    write_final_vcf(deletion_candidates,
                    insertion_candidates,
                    duplication_candidates,
                    breakend_candidates,
                    ref_genome.references,
                    ref_genome.lengths,
                    options)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logging.error(e, exc_info=True)