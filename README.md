# SVPG
[![PyPI version](https://img.shields.io/pypi/v/svpg.svg)](https://pypi.org/project/svpg/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/svpg/badges/version.svg)](https://anaconda.org/bioconda/svpg)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/svpg/badges/license.svg)](https://anaconda.org/bioconda/svpg)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/svpg/badges/platforms.svg)](https://anaconda.org/bioconda/svpg)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/svpg/badges/latest_release_date.svg)](https://anaconda.org/bioconda/svpg)

## Overview
<table style="border-collapse: collapse; border: none; padding: 0; margin: 0; width: 100%;">
  <tr>
    <td style="text-align: center; vertical-align: middle; font-family: monospace; white-space: pre; font-size: 14px; padding: 0; margin: 0;">
<pre style="margin: 0; line-height: 1;">
████ █     █ ████   ████ 
█    █     █ █   █ █     
████  █   █  ████  █ ███ 
   █   █ █   █     █   █ 
████    █    █      ████ 
</pre>
    </td>
    <td vertical-align: middle; padding: 0; margin: 0>
      <div style="margin: 0 auto">
<b>SVPG</b> (Structural Variant detection based on Pangenome Graph) is a computational tool designed for structural variation (SV) detection and efficient pangenome graph augmentation. With the growing availability of long-read sequencing data and pangenome references, SVPG fills a critical gap by enabling accurate SV discovery and scalable integration of new genomes into existing pangenome graphs.
      </div>
    </td>
  </tr>
</table>
<div style="text-align: center; margin-top: 10px;">
  <img src="doc/overview.jpg" alt="SVPG illustration" style="max-width: 100%; height: auto;">
</div>


## Key Features

* **Dual SV detection modes**:

  * **Pangenome-guided mode**:  Extracts SV-supporting reads from BAM files, and realigns a pangenome reference graph. By analyzing the graph alignment's topological and path transition features to detect germline SVs with high precision.
  * **Graph-based mode**: Directly resolves reads-to-graph alignments to discover _de novo_ SVs within haplotype paths of pangenome graph, ideal for conducting reference-bias-free low-frequency/somatic SV discovery without relying on prior SV databases or annotations.
* **High sensitivity and accuracy SV detection**: Demonstrates superior performance in benchmarking against state-of-the-art SV callers across both population-wide germline and individual-specific SVs.
* **Rapid graph augmentation**: Designed to work seamlessly with the graph-call mode, it accelerates pangenome augmentation by nearly an order of magnitude compared to traditional _de novo_ assembly methods on cohorts of dozens of samples, enabling fast and scalable integration of new samples.

## Contents
* [Installation](#installation)
* [Requirements](#requirements)
* [Usage](#usage)
  * [1. Pangenome-Guided SV Detection](#1-pangenome-guided-sv-detection)
  * [2. Graph-Based SV Detection](#2-graph-based-sv-detection)
  * [3. Pangenome Graph Augmentation](#3-pangenome-graph-augmentation)
* [Parameters](#parameters)
* [Limitations](#limitations)
* [Citation](#citation)
* [Contact](#contact)


## Installation

```bash
$ pip install svpg
or
$ conda install svpg
or
$ git clone https://github.com/coopsor/SVPG.git && cd SVPG/ && pip install . 
```

## Requirements
* Python >= 3.10 (tested on v3.10.4)
* pysam >= 0.22 for BAM file processing
* numpy >= 1.26.4 for numerical computing
* scipy >= 1.13.1 for scientific computing
* [pyabpoa](https://github.com/yangao07/abPOA/tree/main/python) >= 1.5.4 for consensus sequence generation

The following tools must be available in your system path (recommend installing via conda):
* [minigraph](https://github.com/lh3/minigraph) >= 0.21 for pangenome graph alignment in pangenome-guided mode
* [mappy](https://github.com/lh3/minimap2/tree/master/python) >= 2.28 for consensus sequence realignment in pangenome-guided mode
* bcftools >= 1.20 for VCFs processing in augmentation mode
* truvari >= 3.1.0 for VCFs merging in augmentation mode

## Usage

### 1. Pangenome-Guided SV Detection
* Pangenome-guided mode requires an input of read-reference alignment results in coordinate-sorted and indexed BAM file. If you start with sequencing reads (e.g., FASTA/FASTQ files), you need to map them to a linear reference genome first.
* SVPG support parallelized and uses 16 threads by default. This value can be adapted using e.g. `-t` 4 as option.
* We evaluated SVPG using the prebuilt human pangenome graph ([v3.1](https://zenodo.org/records/10693675)) constructed from 47 samples. 
* By default, SVPG outputs all SVs supported by more than one read. In pangenome-guided mode, users can according to genotype-assigned variants using `FILTER=PASS` to obtain a more high-confidence SV set.
 In addition, users may manually adjust the minimum read support threshold with the `--min_support`/`-s` parameter based on sequencing depth with the following table for reference. This is particularly useful for ultra-low-coverage datasets (<10×) to preserve recall, as well as for graph-based mode with genotyping is not available.

  | Depth (×) | ONT | HiFi |
  |-----------|-----|------|
  | <10       | 2   | 1    |
  | [10, 20)  | 3   | 2    |
  | [20, 50)  | 4   | 3    |
  | ≥50       | 10  | 4    |

```bash
svpg call --working_dir svpg_out/ --bam sample.bam --ref hg38.fa --gfa pangenome.gfa --read ont
```
The called file `variants.vcf` was saved in the specified working directory. `-o` option can be used to specify the output file name.

### 2. Graph-Based SV Detection
* Graph-based mode requires an input of read-graph alignment results in GAF format. If you start with sequencing reads (e.g., FASTA/FASTQ files), you need to map them to a pangenome. We recommend to produce the alignments using [minigraph]((https://github.com/lh3/minigraph)).
* Since minigraph by default outputs [stable coordinates](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) in [rGFA](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) format, SVPG requires the `--vc` option to be enabled during alignment to support more general GFA formats (e.g., [GraphAligner](https://github.com/maickrau/GraphAligner) alignment result).

```bash
minigraph -cx lr --vc -t 64 pangenome.gfa sample.fasta > sample.gaf 
svpg graph-call --working_dir svpg_out/ --ref hg38.fa --gfa pangenome.gfa --gaf sample.gaf --read ont -s 3
```

* SVPG leverages a pangenome as a panel for filtering germline and population-level SVs, and therefore outputs tumor-only SVs by default. For Tumor/Normal paired analysis, we recommend running the two samples separately and then integrating the results with our script to achieve optimal performance.
```bash
svpg graph-call --working_dir tumor_out/ --ref hg38.fa --gfa pangenome.gfa --gaf tumor.gaf --read hifi -s 3
svpg graph-call --working_dir normal_out/ --ref hg38.fa --gfa pangenome.gfa --gaf normal.gaf --read hifi -s 1
python scripts/vcf_specific.py tumor_out/variants.vcf normal_out/variants.vcf tumor_specific.vcf
```
This procedure selects SVs that are present only in the tumor sample but absent in the matched normal.

### 3. Pangenome Graph Augmentation
SVPG provides a streamlined pipeline to rapidly embed _de novo_ SVs detected from graph-based alignment back into the pangenome graph.
To use this feature, users should place a directory containing the raw sequencing data (e.g., FASTA/FASTQ files) of new samples under the specified `working_dir` path. For example:
```bash
working_dir/
├── sample_1/
│   └── sample_1.fasta
├── sample_2/
│   └── sample_2.fasta
```
SVPG will automatically detect SV in graph-based mode and process these VCFs for graph augmentation, and the output file `augment.gfa` is placed into the given working directory. 
```bash
svpg augment --working_dir svpg_out/ --ref hg38.fa --gfa pangenome.gfa --read hifi
```
Alternatively, you may provide a .tsv file listing the paths to FASTA files of new samples.
For example, the sample.tsv file may look like(sample_1 name ≠ sample_2 name):
`/path/to/sample_1.fasta \n /path/to/sample_2.fasta`
then, run the command `svpg augment --working_dir svpg_out/ --sample_list sample.tsv --ref hg38.fa --gfa pangenome.gfa --read hifi` 

## Parameters
| Parameter               | Description                                                                                                                                                       | Default                                                                            |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------|
| `--working_dir`         | Specify the working directory to store output files.                                                                                                              | Required                                                                           |
| `--bam`                 | Coordinate-sorted and indexed BAM file with aligned long reads.                                                                                                   | Required for `call` mode                                                           |
| `--gaf`                 | GAF file with long reads aligned to the pangenome graph (.gaf).                                                                                                   | Required for `graph-call` mode                                                     |
| `--ref`                 | The reference genome used for pangenome construction (.fa), is also serves as the coordinate system for SVPG’s SV call output.                                    | Required                                                                           |
| `--gfa`                 | Pangenome reference file that the long reads were aligned to (.gfa).                                                                                              | Required                                                                           |
| `--read`                | Type of sequencing reads: `ont` for Oxford Nanopore, `hifi` for PacBio HiFi.                                                                                      | hifi                                                                               |
| `--min_support`/`-s`    | Minimum read support threshold for SV calling. Adjust based on sequencing depth.                                                                                  | 2                                                                                  |
| `--num_threads`/`-t`    | Number of threads to use for parallel processing.                                                                                                                 | 16                                                                                 |
| `--min_mapq`            | Minimum mapping quality for reads to be considered in SV detection.                                                                                               | 20                                                                                 |
| `--min_sv_size`         | Minimum size of SVs to be detected.                                                                                                                               | 50                                                                                 |
| `--max_sv_size`         | Maximum size of SVs to be detected. Set to -1 for unlimited size (recommend for somatic SV of `graph-call` mode).                                                 | 1,000,00                                                                           |
| `--max_merge_threshold` | Maximum distance of SV signals to be merged.                                                                                                                      | 50 for hifi read and 500 for ont read                                              |
| `--ultra_split_size`    | Ignore extremely large BNDs from split alignments unless supported by high enough reads, which may be regarded as false-negative intra-chromosomal translocation. | 1000000                                                                            |
| `--noseq`               | Disable sequence extraction for SVs. Useful for large SVs to save time and disk space.                                                                            | Disabled                                                                           |
| `--types`               | Specify the types of SVs to call: DEL, INS, DUP, INV, BND. Separate multiple types with commas.                                                                   | DEL,INS,DUP,INV,BND                                                                |
| `--contigs`             | Specify the chromosomes list to call SVs (e.g., --contigs chr1 chr2 chrX)'.                                                                                       | All chromosomes                                                                    |   
| `--skip_genotype`       | Skip genotyping step to speed up the process for `call` mode.                                                                                                     | Disabled                                                                           |
| `--realign`             | Realign the noise reads to the reference for more accurate SV sequence inference for `call` mode.                                                                 | Disabled                                                                           |
| `--sample_list`         | Path to a TSV file listing the paths to FASTA files of new samples for `augment` mode.                                                                            | Optional; if not provided, all FASTA files under `working_dir` will be processed.  |
| `--skip_call`           | Skip SV calling step and directly proceed to graph augmentation using existing VCF files in the working directory.                                                | Disabled                                                                           |
| `--out`/`-o`            | Specify the output file name.                                                                                                                                     | `variants.vcf` for `call` and `graph-call` modes, `augment.gfa` for `augment` mode |
| `--version`/`-v`        | Show the version of SVPG.                                                                                                                                         | N/A                                                                                |

## Limitations
* SVPG's pangenome-guided mode relies on minigraph to realign SV signature reads to the pangenome graph. Although this step introduces some overhead, this process is relatively fast: in our tests on the HG002 sample, realignment took approximately 10 minutes for ONT (50×) data and 4 minutes for HiFi (48×) data.
* The `--realign` module provides more accurate breakpoint resolution in graph-hard-alignment regions (for example, [LCRs](https://arxiv.org/html/2509.23057v1#bib.bib20)). On the latest HG002-Q100 benchmark, this module yields measurable performance improvements.
However, it relies on pyabpoa and mappy to perform local re-alignment, which introduces additional computational overhead (e.g., ~1 hour extra for 48× HG002 HiFi data).
As this feature is still experimental, we recommend enabling it in analyses that require base-pair–level breakpoint accuracy.
* The graph-based mode currently does not support genotyping. Users should manually adjust the minimum read support threshold using the `--min_support`/`-s` parameter based on sequencing depth to balance sensitivity and precision.
 
## Citation
Refer to our [paper](https://doi.org/10.1101/2025.07.11.664486) for further details and citation:

Hu, H. et al. SVPG: A pangenome-based structural variant detection approach and rapid augmentation of pangenome graphs with new samples. bioRxiv, 2025.2007.2011.664486 (2025).

## Contact

For questions or support, please open an issue on GitHub or contact the authors at [hhengwork@gmail.com](mailto:hhengwork@gmail.com).
