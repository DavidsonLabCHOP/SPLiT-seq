# SPLiT-seq Demultiplexing Pipeline
This repository contains Davidson Lab's SPLiT-seq demultiplexing pipeline (from FASTQs to gene-per-cell count matrices). All updates are available here.
This tool was created to provide an open source, portable solution for demultiplexing SPLiT-Seq RNA-Seq datasets. 

see https://github.com/paulranum11/SPLiT-Seq_demultiplexing for details on previous versions

for most recent updates and differences between versions see Updates sections

## System Requiremnets

This script has been tested on a linux cluster running Linux CentOS 3.10.0-514.2.2.e17.x86_64 and on a MacBook Pro running macOS High Sierra v10.13.6. NOTE: on macOS systems the options do not work. This problem can be resolved by hardcoding the options you want inside the splitseqdemultiplexing_0.1.4.sh under set default inputs or by installing GNU getopt.

This script is written in bash and python3 and should be portable across a variety of linux systems running the bash shell.

In order to run this software you must install the following dependency packages.

Python3 needs to be installed on your system. Often the executable name of python3 can vary... for example it may appear as python or as python3. This script requires that the executable be python.

Python3 packages: math, os, psutil, argparse, sys, datetime, itertools, re

GNU parallel: https://www.gnu.org/software/parallel/

UMI-tools: https://github.com/CGATOxford/UMI-tools

STAR: https://github.com/alexdobin/STAR in order to run the optional aligment step the STAR rna-seq aligner is required along with an appropriate STAR genome index.

featureCounts: http://subread.sourceforge.net/

Samtools: https://github.com/samtools/samtools

## Getting Started
Download this git repository .zip file or clone this repository using git clone. The downloaded directory will contain three (Round1, Round2, and Round3) barcode files as well as a small example dataset derrived from the 100_CNS_nuclei dataset GEO accession: GSM3017260 (SRR6750041). The full sized datasets can be downloaded from the following European Nucleotide Archive address https://www.ebi.ac.uk/ena/data/view/SRR6750041

The executable file is called splitseqdemultiplex_0.2.3.sh it is written in bash and can be called using bash splitseqdemultiplex_0.2.3.sh (options)

## Options

-n | --numcores # specifies the number of cores you would like to use to parallelize your run.

-v | --version # specifies the version of the demultiplexing utility you would like to run. Input merged to output a single .fastq file with CellIDs and UMIs annotated on each read ID line of the .fastq file. Input split to output a single .fastq file for every single cell sample identified in the input .fastq file with UMIs appended to the read ID line of the .fastq file.

-e | --errors # specifies the number of errors acceptable at each barcode position. The default is set to 1.

-m | --minreads # specifies the minimum number of reads required for a cell to be retained. The default is set to 10.

-1 | --round1barcodes # specifies name of the file containing the barcodes you would like to use for round1. These should be provided as a separate file. See the provided example for formatting reference.

-2 | --round2barcodes # specifies name of the file containing the barcodes you would like to use for round2. These should be provided as a separate file. See the provided example for formatting reference.

-3 | --round3barcodes # specifies name of the file containing the barcodes you would like to use for round3. These should be provided as a separate file. See the provided example for formatting reference.

-f | --fastqF # filepath to the Forward input .fastq file.

-r | --fastqR # filepath to the Reverse input .fastq file.

-o | --outputdir # filepath to the desired output directory.

-t | --targetMemory # define the memory maximum. Processed reads will be saved to memory until this memory maximum is reached. A higher value increases the speed of the script but uses more system memory. Our recommended value is 8000 which equates to 8gb. Higher or lower values will work fine but we suggest using more if your system can support it.

-g | --granularity # the granularity with which you want to get progress updates. Default value is 100000.

-c | --collapseRandomHexamers # when true this option will collapse unique barcode combinations primed with Random Hexamers and OligoDT primers. Because SPLiT-Seq uses both Random Hexamers and OligoDT primers with different barcodes in the same well of the ROUND1 RT step this option is set to true by default.

-a | --align # This is an optional argument. If not included the script will terminate after producing the output .fastq file or files. If star or kallisto are entered as inputs an alignment will be initiated. The star aligner must be used with the -v merged argument and the kallisto aligner may only be used with the -v split option. After alignment, per cell, gene level expression abundance will be calculated and a counts matrix will be produced.

-x | --starGenome # provide the path to the STAR genome index file. Note: STAR options are only relevant if you are doing STAR alignment.

-y | --starGTF # The format of this argument is "GTF /path/to/my/file.gtf" IMPORTANT: The letters GTF must be included before the filepath and the entire argument should be wrapped in quotes as shown in the example. The .gtf file used should be the genome annotation .gtf file that corresponds to the STAR index genome. This is not a file generated by STAR but just a GTF format file downloaded from the same place you downloaded your raw genome. IMPORTANT: Use either a .gtf index using -y or a .saf index using -s. Do not use both.

-s | --geneAnnotationSAF # The format of this argument is "SAF /path/to/my/file.gtf" IMPORTANT: The letters SAF must be included before the filepathand the entire argument should be wrapped in quotes as shown in the example. The .saf file must correspond to the genome used to construct your STAR index. SAF format annotation files are preferred for nuclei sequencing when the user wants to assign both intronic and exonic reads to genes. Prebuilt indexes are provide for Human (GRCh38) and Mouse (GRCm38) from Ensembl.org. IMPORTANT: Use either a .gtf index using -y or a .saf index using -s. Do not use both.

-f | --kallistoIndexIDX # provide the path to your kallisto index file (.idx). Note: Kallisto options are only relevant if you are doing kallisto alignment.

-i | --kallistoIndexFASTA # provide the path to the .fasta file corresponding to your kallisto index file.

Notes: Users may increase the speed of the run by allocating additonal cores using -n and increasing the minimum number of reads required for each cell using -m. Default values for -1 -2 and -3 are the barcodes provided in the splitseq_demultiplexing download: Round1_barcodes_new3.txt, Round2_barcodes_new3.txt and Round3_barcodes_new3.txt. Default values for -f and -r are the provided example .fastq files. The default output directory is results
##Example

# Version 1.0.0+ Update

## Overview
This pipeline provides a robust solution for demultiplexing SPLiT-seq single-cell RNA-seq data and includes optional downstream steps for alignment, feature counting, and UMI processing. By leveraging a C++ implementation for demultiplexing, it ensures high performance while maintaining compatibility with Python for ease of use.
During demultiplexing each "cell" is defined by its unique configuration of SPLiT-Seq round1-3 barcodes.
Optional alignment --align, gene assignment, and counts per gene (per cell) table generation functionality is included.

---

## Description of the Process

1. **Input Splitting**:
   - Splits the input FASTQ files into smaller chunks to enable parallel processing.

2. **Demultiplexing**:
   - Uses a C++ implementation to extract barcodes and UMI information from the reverse reads and transfer them to the forward reads.
   - Filters reads based on a Hamming distance threshold and outputs a merged FASTQ file containing only valid reads.

3. **Performance Metrics (Optional)**:
   - Calculates and reports performance metrics such as the number of valid barcodes and cells passing the minimum read threshold.

4. **Alignment and Post-Processing (Optional)**:
   - Aligns the processed reads to a reference genome using STAR.
   - Counts features using featureCounts.
   - Processes aligned reads for UMI deduplication using UMI-tools.

5. **Final Cleanup**:
   - Removes temporary files and directories to ensure a clean working environment.

---

## Additional prerequisites for versoin 4.0.0+

### Tools and Dependencies
Ensure the following tools are installed and accessible in your environment:

- **Python 3.8+**
- **C++ Compiler** (e.g., `g++` for compiling the demultiplexing library)
- Python Libraries:
  - `argparse`
  - `os`
  - `ctypes`
  - `joblib`

---

## Setup

### Compiling the C++ Demultiplexing Library
1. Place the `barcode_demux.cpp` file in your working directory.
2. Compile the C++ file into a shared library:
   ```bash
   g++ -shared -fPIC -o barcode_demux.so barcode_demux.cpp
   ```
3. Ensure the shared library `barcode_demux.so` is in the same directory as the main Python script.

### Preparing Input Files
1. Ensure the following input files are available:
   - Forward FASTQ file (`inputF.fastq`)
   - Reverse FASTQ file (`inputR.fastq`)
   - Barcode files for Rounds 1, 2, and 3 (e.g., `Round1_barcodes.txt`, etc.)
   - STAR genome directory (if alignment is enabled)
   - Gene annotation file in GTF format (if alignment is enabled)

2. Note the total number of lines in the reverse FASTQ file using:
   ```bash
   wc -l inputR.fastq
   ```

---

## How to Run the Pipeline

### Command
Run the pipeline with the following command:

```bash
python3 splitseqdemultiplex-1.1.0.py \
  -n <num_cores> \
  -e <error_threshold> \
  -m <min_reads_per_cell> \
  -1 <path_to_round1_barcodes> \
  -2 <path_to_round2_barcodes> \
  -3 <path_to_round3_barcodes> \
  -f <forward_fastq_file> \
  -r <reverse_fastq_file> \
  -o <output_directory> \
  -b <reads_per_bin> \
  -l <total_lines_in_reverse_fastq> \
  [-a] \
  [-x <path_to_star_genome>] \
  [-s <path_to_gtf_file>] \
  [-p]
```

### Arguments

- `-n, --numCores`: Number of available CPU cores for splitting files.
- `-e, --errors`: Maximum permissible Hamming distance for barcodes.
- `-m, --minReads`: Minimum number of reads per cell to retain a cell.
- `-1, --round1Barcodes`: Path to the Round 1 barcode file.
- `-2, --round2Barcodes`: Path to the Round 2 barcode file.
- `-3, --round3Barcodes`: Path to the Round 3 barcode file.
- `-f, --fastqF`: Path to the forward FASTQ file.
- `-r, --fastqR`: Path to the reverse FASTQ file.
- `-o, --outputDir`: Directory to store output files.
- `-b, --numReadsBin`: Number of reads processed per bin.
- `-l, --lengthFastq`: Total number of lines in the reverse FASTQ file.
- `-a, --align`: Perform alignment and post-processing (optional).
- `-x, --starGenome`: Path to the STAR genome directory (required if `-a` is used).
- `-s, --geneAnnotation`: Path to the GTF gene annotation file (required if `-a` is used).
- `-p, --performanceMetrics`: Enable performance metrics calculation (optional).

---

## Outputs

### Outputs of demultiplexing segment
- `MergedCells_1.fastq`: Final processed FASTQ file with demultiplexed reads.
- `MergedCells_passing.fastq`: Reads passing the minimum reads-per-cell threshold (if `-p` is enabled).

### Alignment Outputs (if `-a` is enabled, same as before)
- `counts.tsv.gz`: same as previous versions; counts matrix (gene_name x cell_ID).

---

## Example

### Full Pipeline with Alignment
```bash
python3 splitseqdemultiplex-1.1.0.py \
  -n 4 \
  -e 1 \
  -m 10 \
  -1 Round1_barcodes.txt \
  -2 Round2_barcodes.txt \
  -3 Round3_barcodes.txt \
  -f inputF.fastq \
  -r inputR.fastq \
  -o output_directory \
  -b 100000 \
  -l 400000 \
  -a \
  -x /path/to/star/genome \
  -s /path/to/genes.gtf \
  -p
```

---

## Troubleshooting

1. **Demultiplexing Issues**:
   - Check if barcode files are correctly formatted.
   - Ensure input FASTQ files have valid sequences.

2. **Alignment Skipped**:
   - Ensure `-a`, `-x`, and `-s` are provided.

3. **Tool Not Found**:
   - Verify tools are installed and accessible in your `$PATH`.

Previous version (0.#.#) are available by running splitseqdemultiplex-DESIREDVERSION.py **
---

