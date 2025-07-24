import re
import sys
import os
import pickle
import logging
import time
from multiprocessing import Pool
import pysam
from time import strftime, localtime
import subprocess

from svpg.input_parsing import parse_arguments
from svpg.SVCollect import read_bam
from svpg.SVCluster import form_bins, cluster_data
from svpg.SVPan import read_gaf, read_gaf_pan
from svpg.util import read_gfa
from svpg.output_vcf import consolidate_clusters_unilocal, write_final_vcf
from svpg.SVGenotype import genotype
from svpg.graph_augment import augment_pipe

options = parse_arguments()
if options.sub == 'call':
    bam = pysam.AlignmentFile(options.bam, threads=options.num_threads)
ref_genome = pysam.FastaFile(options.ref)

def multi_process(total_len, step, args=None):
    num_threads = int(options.num_threads)
    chunk_size = total_len // num_threads

    analysis_pools = Pool(processes=int(num_threads))
    async_results = []
    for i in range(num_threads):
        start = i * chunk_size
        end = start + chunk_size if i < num_threads - 1 else total_len
        if step == 'read_bam':
            async_results.append(analysis_pools.starmap_async(read_bam, [(args, start, end, options)]))
        elif step == 'read_gaf':
            async_results.append(analysis_pools.starmap_async(read_gaf, [(args[0][start: end], args[1], options)]))
        elif step == 'read_gaf_pan':
            async_results.append(analysis_pools.starmap_async(read_gaf_pan, [(args[0][start: end], args[1], options)]))
        elif step == 'cluster':
            async_results.append(analysis_pools.starmap_async(cluster_data, [(args[0][start:end], args[1])]))
        else:
            async_results.append(analysis_pools.starmap_async(genotype, [(args[0][start:end], args[1], options)]))

    analysis_pools.close()
    analysis_pools.join()
    results = []
    for async_result in async_results:
        result = async_result.get()
        results.extend(result)
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
    for arg in vars(options):
        logging.info("PARAMETER: {0}, VALUE: {1}".format(arg, getattr(options, arg)))

    if options.sub == 'call' or options.sub == 'graph-call':
        if not options.min_support:
            try:
                if options.read == 'ont' or options.read == 'clr':
                    options.min_support = options.depth // 10 + 2
                else:
                    options.min_support = options.depth // 10
            except AttributeError:
                print(
                    "Please specify the min_support or sequencing depth.")
                return

        gfa_node = read_gfa(options.gfa)

        if options.contigs is None:
            options.contigs = [ctg for ctg in ref_genome.references if re.match(r'^(chr)?[0-9XYM]+$', ctg)]

    pan_signatures = []
    if options.sub == 'call':
        logging.info("MODE: call")
        logging.info("INPUT: {0}".format(os.path.abspath(options.bam)))
        logging.info("***************** Collect SV signatures *****************")

        try:
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
        #
        logging.info("****************************** Graph Mapping ******************************")

        # with open(options.working_dir + '/sv_signatures.pkl', 'wb') as temp:
        #     pickle.dump(bam_signatures, temp)
        # with open(options.working_dir +'/sv_signatures.pkl', 'rb') as f:
        #     bam_signatures = pickle.load(f)
        if not os.path.exists(options.working_dir + '/signatures.fa'):
            fasta_file = open(options.working_dir + '/signatures.fa', 'a')
            signature_bin, bin_depth = form_bins(bam_signatures, 1000)
            signature_bin_remap = [cluster for cluster in signature_bin if len(cluster) >= options.min_support]
            for cluster in signature_bin_remap:
                for sig in cluster:
                    pos_ref = str(sig.contig) + ':' + str(sig.start) + ':' + str(sig.end)
                    # adjac_distance = max(min(5000, sig.svlen*3), 2000)
                    adjac_distance = 2000
                    if sig.signature == 'suppl':
                        read_seq = sig.read_seq
                    elif sig.signature == 'merged_indel':
                        read_seq = sig.read_seq
                    else:
                        if sig.pos_read < adjac_distance:
                            read_seq = sig.read_seq[0:sig.pos_read + sig.svlen + adjac_distance]
                        else:
                            read_seq = sig.read_seq[sig.pos_read - adjac_distance:sig.pos_read + sig.svlen + adjac_distance] if sig.pos_read else sig.read_seq

                    ref_suppl1, ref_suppl2 = '', ''
                    if sig.pos_read < adjac_distance:
                        try:
                            ref_suppl1 = ref_genome.fetch(sig.contig, sig.start - adjac_distance, sig.start - sig.pos_read)
                        except ValueError:
                            ref_suppl1 = ''
                        sig.signature = 'short'
                    if sig.type == 'DEL' and sig.pos_read + adjac_distance > len(sig.read_seq):    # 窗口右端超出read末端
                        try:
                            ref_suppl2 = ref_genome.fetch(sig.contig, sig.end + (len(sig.read_seq) - sig.pos_read),
                                                          sig.end + adjac_distance)
                        except ValueError:
                            ref_suppl2 = ref_genome.fetch(sig.contig, sig.end + (len(sig.read_seq) - sig.pos_read),
                                                          len(ref_genome.fetch(sig.contig)))
                        sig.signature = 'short'
                    elif sig.type == 'INS' and sig.pos_read + sig.svlen + adjac_distance > len(sig.read_seq):
                        try:
                            ref_suppl2 = ref_genome.fetch(sig.contig, sig.start + len(sig.read_seq) - sig.pos_read - sig.svlen,
                                                      sig.start + adjac_distance)
                        except ValueError:
                            ref_suppl2 = ref_genome.fetch(sig.contig, sig.start + len(sig.read_seq) - sig.pos_read - sig.svlen,
                                                          len(ref_genome.fetch(sig.contig)))
                        sig.signature = 'short'
                    read_seq = ref_suppl1 + read_seq + ref_suppl2
                    fasta_file.write(f'>{sig.read_name}@{sig.type}@{pos_ref}@{sig.signature}\n{read_seq}\n')
            fasta_file.close()

            if options.read == 'hifi':
                os.system(
                        'minigraph -t 128 -cx asm --vc --secondary yes  {0} {1}/signatures.fa > {1}/signatures.gaf'.format(options.gfa, options.working_dir))
            else:
                os.system(
                        'minigraph -t 128 -cx lr --vc --secondary yes {0} {1}/signatures.fa > {1}/signatures.gaf'.format(options.gfa, options.working_dir))

        logging.info("*************** Collect signatures from pangenome-reference ***************")

        with open(options.working_dir + '/signatures.gaf', 'rb') as f:
            for chunk_index, lines in enumerate(read_in_chunks(f, chunk_size=200000000)):
                logging.info(f"Processing chunks {chunk_index + 1}")
                pan_signatures.extend(multi_process(len(lines), 'read_gaf', (lines, gfa_node)))
                logging.info(f"Processed chunks {chunk_index + 1}")

    elif options.sub == 'graph-call':
        logging.info("MODE: graph-call")
        logging.info("INPUT: {0}".format(os.path.abspath(options.gaf)))
        logging.info("*************** Collect SV signatures from pangenome ***************")

        with open(options.gaf, 'rb') as f:
            for chunk_index, lines in enumerate(read_in_chunks(f, chunk_size=200000000)):
                logging.info(f"Processing chunk {chunk_index+1}")
                pan_signatures.extend(multi_process(len(lines), 'read_gaf_pan', (lines, gfa_node)))
                logging.info(f"Processed chunks {chunk_index+1}")

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
                    if entry.is_dir() and entry.name.startswith("HG"):
                        # Assuming FASTA file is directly inside the sample directory and named as {prefix}.fasta
                        fasta_path_in_dir = os.path.join(entry.path, f"{entry.name}.fasta")
                        if os.path.exists(fasta_path_in_dir):
                            sample_paths_to_process.append(fasta_path_in_dir)
                        else:
                            logging.warning(
                                f"Expected FASTA file {fasta_path_in_dir} not found. Skipping {entry.name}.")
                            raise RuntimeError(f"FASTA file not found for sample: {entry.name}")

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

                        logging.info(f"Start call SVs from sample {prefix}")
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
                                cmd_align = f"minigraph -t128 -cxasm --vc --secondary yes {options.gfa} {fasta_file_name} > {gaf_file}"
                            else:
                                cmd_align = f"minigraph -t128 -cxlr --vc --secondary yes {options.gfa} {fasta_file_name} > {gaf_file}"
                            subprocess.run(cmd_align, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                           text=True)

                        var_file = options.vcf_out
                        cmd_call = [
                            "python", __file__, "graph-call",
                            "--read", options.read,
                            "-s", str(support),
                            "--working_dir", './',  # This should refer to the current sample directory
                            "--ref", options.ref,
                            "--gfa", options.gfa,
                            "--gaf", gaf_file,
                            "-o", var_file,
                            "--raw_fasta", fasta_file_name,
                            "--min_sv_size", str(options.min_sv_size),
                            "--max_sv_size", str(options.max_sv_size)
                        ]
                        try:
                            subprocess.run(cmd_call, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
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

    # with open(options.working_dir + '/pan_signatures.pkl', 'wb') as temp:
    #     pickle.dump(pan_signatures, temp)
    # with open(options.working_dir+'/pan_signatures.pkl', 'rb') as f:
    #     pan_signatures = pickle.load(f)
    deletion_signatures = [ev for ev in pan_signatures if ev.type == "DEL"]
    insertion_signatures = [ev for ev in pan_signatures if ev.type == "INS"]
    logging.info("Found {0} signatures for deleted regions.".format(len(deletion_signatures)))
    logging.info("Found {0} signatures for inserted regions.".format(len(insertion_signatures)))

    logging.info("**************************** Cluster signatures ***************************")

    filted_signature = []
    bin_depth_list = []
    for element_signature in [deletion_signatures, insertion_signatures]:
        if not element_signature:
            continue
        signature_bin, bin_depth = form_bins(element_signature, 1000)
        bin_depth_list.append(bin_depth)
        signature_clusters = []
        signature_clusters.extend(multi_process(len(signature_bin), 'cluster', (signature_bin, bin_depth)))
        entropies, node_ls_collection = [], []
        for cluster in signature_clusters:
            entropie_flag = 0
            for current, next_item in zip(cluster[:-1], cluster[1:]):
                if current.node_ls and next_item.node_ls:
                    if set(current.node_ls) != set(next_item.node_ls):
                        entropie_flag = 1
                        break
                else:
                    if current.svlen != next_item.svlen or current.start != next_item.start:
                        entropie_flag = 1
                        break
            entropies.append(entropie_flag)
        for ne, cluster in zip(entropies, signature_clusters):
            if len(cluster) >= options.min_support:
                filted_signature.append(cluster)
            else:
                if options.min_support > 2:
                    if options.min_support - 1 == len(cluster) and ne == 0:
                        filted_signature.append(cluster)

        logging.info(f"Generated {len(signature_clusters)} signature clusters")

    chrom_results = {}
    for sig in filted_signature:
        chrom_results.setdefault(sig[0].contig, []).append(sig)

    logging.info("********************************** SVCALL *********************************")

    sv_candidate = []
    fasta_file = pysam.FastaFile(options.raw_fasta) if getattr(options, "raw_fasta", None) else None

    for contig in chrom_results:
        if contig in options.contigs:
            ref_chrom_seq = ref_genome.fetch(contig)
            sv_candidate.extend(sorted(
                consolidate_clusters_unilocal(chrom_results[contig], ref_chrom_seq, options, fasta_file),
                key=lambda cluster: (cluster.contig, cluster.start)))

    deletion_candidates = [i for i in sv_candidate if i.type == 'DEL']
    insertion_candidates = [i for i in sv_candidate if i.type == 'INS']
    logging.info("Final deletion candidates: {0}".format(len(deletion_candidates)))
    logging.info("Final insertion candidates: {0}".format(len(insertion_candidates)))

    if options.sub == 'call' and not options.skip_genotype:
        logging.info("********************************* GENOTYPE ********************************")
        logging.info("Genotyping deletions..")
        deletion_candidates = multi_process(len(deletion_candidates), 'genotype', (deletion_candidates, "DEL"))
        logging.info("Genotyping insertions..")
        insertion_candidates = multi_process(len(insertion_candidates), 'genotype', (insertion_candidates, "INS"))

    types_to_output = [entry.strip() for entry in options.types.split(",")]
    write_final_vcf(deletion_candidates,
                    insertion_candidates,
                    ref_genome.references,
                    ref_genome.lengths,
                    types_to_output,
                    options.working_dir,
                    options.out)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logging.error(e, exc_info=True)