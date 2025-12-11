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
                        version='svpg v1.3.0')

    subparsers = parser.add_subparsers(help='mode', dest='sub')
    parser.set_defaults(sub='call')

    parser_bam = subparsers.add_parser('call',
                                       help='Pangenome-guided SV detection')
    parser_bam.add_argument('--working_dir',
                            type=os.path.abspath,
                            help='Specify the working directory to store output files.' )
    parser_bam.add_argument('--bam',
                            type=str,
                            help='Coordinate-sorted and indexed BAM file with aligned long reads.')
    parser_bam.add_argument('--ref',
                            type=str,
                            help='The reference genome used for pangenome construction (.fa), is also serves as the coordinate system for SVPG’s SV call output.')
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
                            choices=['hifi', 'ont'],
                            default='hifi',
                            help="Type of sequencing reads: `ont` for Oxford Nanopore, `hifi` for PacBio HiFi. ")
    parser_bam.add_argument('--min_mapq',
                            type=int,
                            default=20,
                            help='Minimum mapping quality for reads to be considered in SV detection.')
    parser_bam.add_argument('--min_sv_size',
                            type=int,
                            default=50,
                            help='Minimum size of SVs to be detected.')
    parser_bam.add_argument('--max_sv_size',
                            type=int,
                            default=100000,
                            help='Maximum SV size to detect include sequence information. Set to -1 for unlimited size.')
    parser_bam.add_argument('--max_merge_threshold',
                            type=int,
                            default=None,
                            help='Maximum distance of SV signals to be merged.')
    parser_bam.add_argument('--ultra_split_size',
                            type=int,
                            default=1000000,
                            help='Ignore extremely large BNDs from split alignments unless supported by high enough reads,\
                                  which may be regarded as false-negative intra-chromosomal translocation')
    parser_bam.add_argument('--noseq',
                            action='store_true',
                            help='Disable sequence extraction for SVs. Useful for large SVs to save time and disk space.')
    parser_bam.add_argument('--realign',
                            action='store_true',
                            help='Realign the noise reads to the reference for more accurate SV sequence inference')
    parser_bam.add_argument('-s', '--min_support',
                            type=int,
                            default=2,
                            help='Minimum read support threshold for SV calling. Adjust based on sequencing depth.')
    parser_bam.add_argument('--types',
                            type=str,
                            default="DEL,INS,DUP,INV,BND",
                            help='SV types to include in output VCF (default: %(default)s). \
                                  Give a comma-separated list of SV types, like "DEL,INS"')
    parser_bam.add_argument('--contigs',
                            type=str,
                            nargs='*',
                            help='Specify the chromosomes list to call SVs (e.g., --contigs chr1 chr2 chrX)')
    parser_bam.add_argument('--skip_genotype',
                            action='store_true',
                            help='Skip genotyping step to speed up the processing.')

    ##########################################################
    parser_gaf = subparsers.add_parser('graph-call',
                                       help='Pangenome-based de novo SV detection')
    parser_gaf.add_argument('--working_dir',
                            type=os.path.abspath,
                            help='Specify the working directory to store output files.')
    parser_gaf.add_argument('--ref',
                            type=str,
                            help='The reference genome used for pangenome construction (.fa), is also serves as the coordinate system for SVPG’s SV call output.')
    parser_gaf.add_argument('--gfa',
                            type=str,
                            help='Pangenome reference file that the long reads were aligned to (.gfa)')
    parser_gaf.add_argument('--gaf',
                            type=str,
                            help='GAF file that aligns to the pangenome reference (.gaf)')
    parser_gaf.add_argument('-o', '--out',
                            type=str,
                            default='variants.vcf',
                            help='Specify the output file name.')
    parser_gaf.add_argument('-t', '--num_threads',
                            type=int,
                            default=16,
                            help='Number of threads to use for parallel processing.')
    parser_gaf.add_argument('--read',
                            type=str,
                            default='hifi',
                            help='Type of sequencing reads: `ont` for Oxford Nanopore, `hifi` for PacBio HiFi.')
    parser_gaf.add_argument('--min_mapq',
                            type=int,
                            default=20,
                            help='Minimum mapping quality for reads to be considered in SV detection.')
    parser_gaf.add_argument('--max_merge_threshold',
                            type=int,
                            default=500,
                            help='Maximum distance of SV signals to be merged')
    parser_gaf.add_argument('--min_sv_size',
                            type=int,
                            default=50,
                            help='Minimum size of SVs to be detected.')
    parser_gaf.add_argument('--max_sv_size',
                            type=int,
                            default=100000,
                            help='Maximum size of SVs to be detected. Set to -1 for unlimited size (recommend somatic SV).')
    parser_gaf.add_argument('--ultra_split_size',
                            type=int,
                            default=1000000,
                            help='Ignore extremely large BNDs from split alignments unless supported by high enough reads,\
                                  which may be regarded as false-negative intra-chromosomal translocation')
    parser_gaf.add_argument('--noseq',
                            action='store_true',
                            help='Disable sequence extraction for SVs. Useful for large SVs to save time and disk space.')
    parser_gaf.add_argument('-s', '--min_support',
                            type=int,
                            default=2,
                            help='Minimum read support threshold for SV calling. Adjust based on sequencing depth.')
    parser_gaf.add_argument('--types',
                            type=str,
                            default="DEL,INS,DUP,INV,BND",
                            help='SV types to include in output VCF (default: %(default)s). \
                                  Give a comma-separated list of SV types, like "DEL,INS"')
    parser_gaf.add_argument('--contigs',
                            type=str,
                            nargs='*',
                            help='Specify the chromosomes list to call SVs (e.g., --contigs chr1 chr2 chrX)')

    ##########################################################
    parser_augment = subparsers.add_parser('augment',
                                           help='Pangenome graph augmentation pipeline')
    parser_augment.add_argument('--working_dir',
                                type=os.path.abspath,
                                help='Specify the working directory to store output files.')
    parser_augment.add_argument('--ref',
                                type=str,
                                help='The reference genome used for pangenome construction (.fa), is also serves as the coordinate system for SVPG’s SV call output.')
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
                                help='Number of threads to use for parallel processing.')
    parser_augment.add_argument('--read',
                                type=str,
                                default='hifi',
                                help='Type of sequencing reads: `ont` for Oxford Nanopore, `hifi` for PacBio HiFi.')
    parser_augment.add_argument('--sample_list',
                                type=str,
                                default='',
                                help='Path to a TSV file listing the paths to FASTA files of new samples. if not provided, all FASTA files under `working_dir` will be processed.\
                                For example, the sample.tsv file may look like(sample_1 name ≠ sample_2 name): /path/to/sample_1.fasta\n/path/to/sample_2.fasta')
    parser_augment.add_argument('--min_mapq',
                                type=int,
                                default=20,
                                help='Minimum mapping quality for reads to be considered in SV detection.')
    parser_augment.add_argument('--min_sv_size',
                                type=int,
                                default=50,
                                help='Minimum size of SVs to be detected.')
    parser_augment.add_argument('--max_sv_size',
                                type=int,
                                default=-1,
                                help='Minimum size of SVs to be detected. Set to -1 for unlimited size (recommend somatic SV).')
    parser_augment.add_argument('--skip_call',
                                action='store_true',
                                help='Skip SV calling step and directly proceed to graph augmentation using existing VCF files in the working directory. ')

    return parser.parse_args(arguments)
