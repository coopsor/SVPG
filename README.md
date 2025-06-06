# SVPG

## Overview

SVPG (Structural Variant detection based on Pangenome Graph) is a computational tool designed for structural variation (SV) detection and efficient pangenome graph augmentation. With the growing availability of long-read sequencing data and pangenome references, SVPG fills a critical gap by enabling accurate SV discovery and scalable integration of new genomes into existing pangenome graphs.

## Key Features

* **Dual SV Detection Modes**:

  * **Pangenome-Guided Mode**:  Extracts SV-supporting reads from BAM files, converts them into signature sequences, and realigns them against a pangenome reference graph. By analyzing the graph alignment's topological and path transition features to detect SVs with high precision.
  * **Graph-Based Mode**: Directly resolves read-to-graph alignments to discover _de novo_ SVs within haplotype paths of pangenome graph, ideal for conducting reference-bias-free low-frequency SV discovery without relying on prior SV databases or annotations.
* **High Sensitivity and Accuracy**: Demonstrates superior performance in benchmarking against state-of-the-art SV callers across both population-wide germline and individual-specific SVs.
* **Rapid Graph Augmentation**: Designed to work seamlessly with the graph-call mode, it accelerates pangenome augmentation by nearly an order of magnitude compared to traditional _de novo_ assembly methods on cohorts of dozens of samples, enabling fast and scalable integration of new samples.
## Installation

```bash
$ git clone https://github.com/coopsor/SVPG.git && cd SVPG/ && pip install . 
```

## Requirements
* Python >= 3.10.4 (tested on v3.10.4)
* pysam >= 0.22
* numpy >= 1.26.4 

The following tools must be available in your system path (recommend installing via conda):
* minigraph >= 0.21
* bcftools >= 1.20
* truvari >= 3.1.0

## Usage

### 1. Pangenome-Guided SV Detection
* SVPG support parallelized and uses 16 threads by default. This value can be adapted using e.g. `-t` 4 as option.
* We evaluated SVPG using the prebuilt human pangenome graph ([v3.1](https://zenodo.org/records/10693675)) constructed from 47 samples, which represented the latest version at the time of testing prior to public release. 
* The following table provides recommended `--min_support` values to filter out low-quality SVs under different sequencing depths for **ONT** and **HiFi** platforms. Alternatively, users can specify the sequencing depth with `-d`(`--depth`), and SVPG will automatically estimate an appropriate minimum support threshold.

  | Depth Range  | ONT | HiFi |
  |--------------|-----|------|
  | 5 to <10     | 2   | 1    |
  | 10 to <20    | 3   | 2    |
  | 20 to <50    | 4   | 3    |
  | ≥50          | 10  | 4    |

```bash
svpg call --working_dir svpg_out --bam sample.bam --ref hg38.fa --gfa pangenome.gfa --read ont -s min_support
```

### 2. Graph-Based SV Detection
* Graph-based mode requires an input of read-graph alignment results in GAF format. If you start with sequencing reads (Fasta or Fastq format), you may use minigraph to map them to a pangenome reference.
* Since minigraph by default outputs [stable coordinates](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) in [rGFA](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) format, SVPG requires the `--vc` option to be enabled during alignment to support more general GFA formats (e.g., [GraphAligner](https://github.com/maickrau/GraphAligner) alignment result)

```bash
svpg graph-call --working_dir svpg_out --ref hg38.fa --gfa pangenome.gfa --gaf sample.gaf --read ont -s min_support
```

### 3. Pangenome Graph Augmentation
SVPG provides a streamlined pipeline to rapidly embed _de novo_ SVs detected from graph-based alignment back into the pangenome graph.
To use this feature, users should place a directory containing the raw sequencing data (e.g., FASTQ files) of new samples under the specified `working_dir` path. For example:
```bash
working_dir/
├── sample_1/
│   └── sample_1.fasta
├── sample_2/
│   └── sample_2.fasta
```
SVPG will automatically detect SV and process these files for graph augmentation:
```bash
svpg augment --working_dir svpg_out --ref hg38.fa --gfa pangenome.gfa --read hifi
```
Alternatively, you may provide a .tsv file listing the paths to FASTA files of new samples.
For example, the sample.tsv file may look like(sample_1 name ≠ sample_2 name):
`/path/to/sample_1.fasta \n /path/to/sample_2.fasta`
then, run the command `svpg augment --working_dir svpg_out/ --sample_list sample.tsv --ref hg38.fa --gfa pangenome.gfa --read hifi` 

## Limitations
* SVPG's pangenome-guided mode relies on minigraph to realign SV signature reads to the reference pangenome graph. Although this step introduces some overhead, this process is relatively fast: in our tests on the HG002 sample, realignment took approximately 10 minutes for ONT (50×) data and 4 minutes for HiFi (48×) data.
* SVPG is not a dedicated somatic SV caller like [Severus](https://github.com/KolmogorovLab/Severus) or [Savana](https://github.com/cortes-ciriano-lab/savana), and therefore may have limited ability to detect complex BND events. However, SVPG enable filtering out common germline SVs using the pangenome significantly enhances the detection of somatic-specific indels. This finding highlights a promising direction for somatic SV research. It is possible to construct personalized or population-specific pangenome references from matched normal (adjacent) tissues could enable more accurate detection of somatic SV.
* The graph-based SV calling mode currently does not support genotyping. However, genotypes can be inferred based on the number of reads supporting the SV and the reference paths in the graph, the functionality will be added in future versions.

## Citation


## Contact

For questions or support, please open an issue on GitHub or contact the authors at [hhengwork@gmail.com](mailto:hhengwork@gmail.com).
