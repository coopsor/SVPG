import sys
import os
import argparse


def parse_arguments(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVPG - Structural variant detection based on pangenome graph""")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('-v', '--version',
                        action='version',
                        version='svpg v1.2.0')

    subparsers = parser.add_subparsers(help='mode', dest='sub')
    parser.set_defaults(sub='call')

    parser_bam = subparsers.add_parser('call',
                                       help='Pangenome-guided SV detection')
    parser_bam.add_argument('--working_dir',
                            type=os.path.abspath,
                            help='Working and output directory. \
                                  Existing files in the directory are overwritten. \
                                  If the directory does not exist, it is created.')
    parser_bam.add_argument('--bam',
                            type=str,
                            help='Coordinate-sorted and indexed BAM file with aligned long reads')
    parser_bam.add_argument('--ref',
                            type=str,
                            help='Liner reference genome file that the long reads were aligned to (.fa)')
    parser_bam.add_argument('-o', '--out',
                            type=str,
                            default='variants.vcf',
                            help='VCF output file name')
    parser_bam.add_argument('--gfa',
                            type=str,
                            help='Pangenome reference file that the long reads were aligned to (.gfa)')
    parser_bam.add_argument('-t', '--num_threads',
                            type=int,
                            default=16,
                            help='Number of threads to use')
    parser_bam.add_argument('--read',
                            type=str,
                            default='ont',
                            help='Platform type for sequencing data, either "ont", "hifi" or "clr" (default: %(default)s)')
    parser_bam.add_argument('--min_mapq',
                            type=int,
                            default=20,
                            help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                  Reads with a lower mapping quality are ignored.')
    parser_bam.add_argument('--min_sv_size',
                            type=int,
                            default=40,
                            help='Minimum SV size to detect (default: %(default)s). \
                                  SVPG can potentially detect events of any size but is limited by the \
                                  signal-to-noise ratio in the input alignments. That means that more \
                                  accurate reads and alignments enable the detection of smaller events. \
                                  For current PacBio or Nanopore data, we would recommend a minimum size \
                                  of min_sv_sizebp or larger.')
    parser_bam.add_argument('--max_sv_size',
                            type=int,
                            default=-1,
                            help='Maximum SV size to detect include sequence information. Set to -1 for unlimited size.')
    parser_bam.add_argument('--noseq',
                            action='store_true',
                            default=False,
                            help='Omit sequence information in VCF. Recommended for somatic SV analysis with many large SVs.')
    parser_bam.add_argument('-d', '--depth',
                            type=int,
                            help='Sequencing depth of this dataset')
    parser_bam.add_argument('-s', '--min_support',
                            type=int,
                            help='Minimal number of supporting reads for one SV event')
    parser_bam.add_argument('--types',
                            type=str,
                            default="DEL,INS",
                            help='SV types to include in output VCF (default: %(default)s). \
                                  Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                  INS (novel insertions)')
    parser_bam.add_argument('--contigs',
                            type=str,
                            nargs='*',
                            help='Specify the chromosomes list to call SVs (e.g., --contigs chr1 chr2 chrX)')
    parser_bam.add_argument('--skip_genotype',
                            action='store_true',
                            help='Skip genotyping and only call SVs')

    ##########################################################
    parser_gaf = subparsers.add_parser('graph-call',
                                       help='Pangenome-based de novo SV detection')
    parser_gaf.add_argument('--working_dir',
                            type=os.path.abspath,
                            help='Working and output directory. \
                                  Existing files in the directory are overwritten. \
                                  If the directory does not exist, it is created.')
    parser_gaf.add_argument('--ref',
                            type=str,
                            help='Liner reference genome file that the long reads were aligned to (.fa)')
    parser_gaf.add_argument('--gfa',
                            type=str,
                            help='Pangenome reference file that the long reads were aligned to (.gfa)')
    parser_gaf.add_argument('--gaf',
                            type=str,
                            help='GAF file that aligns to the pangenome reference (.gaf)')
    parser_gaf.add_argument('-o', '--out',
                            type=str,
                            default='variants.vcf',
                            help='VCF output file name')
    parser_gaf.add_argument('-t', '--num_threads',
                            type=int,
                            default=16,
                            help='Number of threads to use')
    parser_gaf.add_argument('--read',
                            type=str,
                            default='ont',
                            help='Platform type for sequencing data, either "ont", "hifi" or "clr"  (default: %(default)s)')
    parser_gaf.add_argument('--raw_fasta',
                            type=str,
                            help='Raw fasta file that the long reads. Since the GAF format does not store \
                                  complete sequences information, providing the original reads\
                                  sequencing data through the --raw_fasta parameter becomes essential \
                                  in graph augmentation mode or when forced output of ALT information \
                                  for split alignments is required. ')
    parser_gaf.add_argument('--min_mapq',
                            type=int,
                            default=20,
                            help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                  Reads with a lower mapping quality are ignored.')
    parser_gaf.add_argument('--min_sv_size',
                            type=int,
                            default=40,
                            help='Minimum SV size to detect (default: %(default)s). \
                                  SVPG can potentially detect events of any size but is limited by the \
                                  signal-to-noise ratio in the input alignments. That means that more \
                                  accurate reads and alignments enable the detection of smaller events. \
                                  For current PacBio or Nanopore data, we would recommend a minimum size \
                                  of 40bp or larger.')
    parser_gaf.add_argument('--max_sv_size',
                            type=int,
                            default=-1,
                            help='Maximum SV size to detect include sequence information. Set to -1 for unlimited size.')
    parser_gaf.add_argument('--noseq',
                            action='store_true',
                            default=False,
                            help='Omit sequence information in VCF. Recommended for somatic SV analysis with many large SVs.')
    parser_gaf.add_argument('-d', '--depth',
                            type=int,
                            help='Sequencing depth of this dataset')
    parser_gaf.add_argument('-s', '--min_support',
                            type=int,
                            help='Minimal number of supporting reads for one SV event')
    parser_gaf.add_argument('--types',
                            type=str,
                            default="DEL,INS",
                            help='SV types to include in output VCF (default: %(default)s). \
                                  Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                  INS (novel insertions)')
    parser_gaf.add_argument('--contigs',
                            type=str,
                            nargs='*',
                            help='Specify the chromosomes list to call SVs (e.g., --contigs chr1 chr2 chrX)')

    ##########################################################
    parser_augment = subparsers.add_parser('augment',
                                           help='Pangenome graph augmentation pipeline')
    parser_augment.add_argument('--working_dir',
                                type=os.path.abspath,
                                help='Working and output directory. \
                                      Existing files in the directory are overwritten. \
                                      If the directory does not exist, it is created.')
    parser_augment.add_argument('--ref',
                                type=str,
                                help='Liner reference genome file that the long reads were aligned to (.fa)')
    parser_augment.add_argument('--gfa',
                                type=str,
                                help='Pangenome reference file that the long reads were aligned to (.gfa)')
    parser_augment.add_argument('-o', '--out',
                                type=str,
                                default='augment.gfa',
                                help='Augmented GFA output file name')
    parser_augment.add_argument('--vcf_out',
                                type=str,
                                default='variants.vcf',
                                help='VCF output file name')
    parser_augment.add_argument('-t', '--num_threads',
                                type=int,
                                default=16,
                                help='Number of threads to use')
    parser_augment.add_argument('--read',
                                type=str,
                                default='ont',
                                help='Platform type for sequencing data, either "ont", "hifi" or "clr"  (default: %(default)s)')
    parser_augment.add_argument('--sample_list',
                                type=str,
                                default='',
                                help='Provide a .tsv file listing the paths to FASTA files of new samples. \
                                For example, the sample.tsv file may look like(sample_1 name ≠ sample_2 name): /path/to/sample_1.fasta\n/path/to/sample_2.fasta')
    parser_augment.add_argument('--min_mapq',
                                type=int,
                                default=20,
                                help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                        Reads with a lower mapping quality are ignored.')
    parser_augment.add_argument('--min_sv_size',
                                type=int,
                                default=40,
                                help='Minimum SV size to detect (default: %(default)s). \
                                      SVPG can potentially detect events of any size but is limited by the \
                                      signal-to-noise ratio in the input alignments. That means that more \
                                      accurate reads and alignments enable the detection of smaller events. \
                                      For current PacBio or Nanopore data, we would recommend a minimum size \
                                      of 40bp or larger.')
    parser_augment.add_argument('--max_sv_size',
                                type=int,
                                default=-1,
                                help='Maximum SV size to detect include sequence information. Set to -1 for unlimited size.')
    parser_augment.add_argument('--skip_call',
                                action='store_true',
                                help='Skip call SVs and only graph augment')

    return parser.parse_args(arguments)
