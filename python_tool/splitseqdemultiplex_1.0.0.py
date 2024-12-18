#!/usr/bin/env python
# Updated SPLiT-seq script with C++ demultiplexing and full pipeline integration

import Splitseq_fun_lib
import os
import argparse
import ctypes

# Load the C++ shared library
lib_path = "./barcode_demux.so"
barcode_demux = ctypes.CDLL(lib_path)

# Define the argument types and return type for the C++ function
barcode_demux.demultiplex.argtypes = [
    ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
    ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
    ctypes.c_int, ctypes.c_int
]
barcode_demux.demultiplex.restype = None

# Wrapper for the C++ demultiplex function
def run_demultiplex(inputF, inputR, output, round1, round2, round3, error_threshold, num_reads_bin):
    barcode_demux.demultiplex(
        inputF.encode('utf-8'), inputR.encode('utf-8'), output.encode('utf-8'),
        round1.encode('utf-8'), round2.encode('utf-8'), round3.encode('utf-8'),
        error_threshold, num_reads_bin
    )

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--numCores', required=True, help='Number of available CPUs for splitting files.')
parser.add_argument('-e', '--errors', required=True, help='Max number of errors permissible per barcode.')
parser.add_argument('-m', '--minReads', required=True, help='Minimum number of reads per cell for retention.')
parser.add_argument('-1', '--round1Barcodes', required=True, help='Path to the Round 1 barcode file.')
parser.add_argument('-2', '--round2Barcodes', required=True, help='Path to the Round 2 barcode file.')
parser.add_argument('-3', '--round3Barcodes', required=True, help='Path to the Round 3 barcode file.')
parser.add_argument('-f', '--fastqF', required=True, help='Path to the forward FASTQ file.')
parser.add_argument('-r', '--fastqR', required=True, help='Path to the reverse FASTQ file.')
parser.add_argument('-o', '--outputDir', required=True, help='Directory to store output files.')
parser.add_argument('-b', '--numReadsBin', required=True, help='Number of reads processed per bin.')
parser.add_argument('-l', '--lengthFastq', required=True, help='Total number of lines in the reverse FASTQ file.')
parser.add_argument('-a', '--align', required=False, action='store_true', help='Perform alignment step.')
parser.add_argument('-x', '--starGenome', required=False, help='Path to the STAR genome directory.')
parser.add_argument('-s', '--geneAnnotation', required=False, help='Path to the gene annotation file (GTF format).')
parser.add_argument('-p', '--performanceMetrics', required=False, action='store_true', help='Enable performance metrics.')
args = parser.parse_args()

#####################################
# STEP 1: Split Input FASTQ Files   #
#####################################
print("Step 1: Splitting input FASTQ files for sequential processing...")

# Create the output directory if it doesn't exist
if not os.path.exists(args.outputDir):
    os.makedirs(args.outputDir)

# Split forward and reverse FASTQ files
Splitseq_fun_lib.split_fastqF_fun(int(args.numCores), args.fastqF, int(args.lengthFastq))
Splitseq_fun_lib.split_fastqR_fun(int(args.numCores), args.fastqR, int(args.lengthFastq))

##############################################
# STEP 2: Sequential Demultiplexing with C++ #
##############################################
print("Step 2: Demultiplexing FASTQ files sequentially using C++ library...")

for i in range(int(args.numCores)):
    run_demultiplex(
        inputF=f"split_fastq_F_{i}",
        inputR=f"split_fastq_R_{i}",
        output=f"{args.outputDir}/demux_{i}.fastq",
        round1=args.round1Barcodes,
        round2=args.round2Barcodes,
        round3=args.round3Barcodes,
        error_threshold=int(args.errors),
        num_reads_bin=int(args.numReadsBin)
    )

##############################################
# STEP 3: Merge Demultiplexed Outputs        #
##############################################
print("Step 3: Merging demultiplexed FASTQ files...")

merged_file = os.path.join(args.outputDir, "MergedCells_1.fastq")
with open(merged_file, "w") as outfile:
    for i in range(int(args.numCores)):
        split_file = f"{args.outputDir}/demux_{i}.fastq"
        with open(split_file, "r") as infile:
            outfile.write(infile.read())
        os.remove(split_file)

##############################################
# STEP 4: Extract Performance Metrics        #
##############################################
if args.performanceMetrics:
    print("Step 4: Extracting performance metrics...")
    import DemultiplexUsingBarcodes_New_V2
    DemultiplexUsingBarcodes_New_V2.calc_demux_results(
        outputDir=args.outputDir,
        performanceMetrics=True,
        readsPerCellThreshold=int(args.minReads)
    )
else:
    print("Step 4: Skipping performance metrics extraction.")

##############################################
# STEP 5: Alignment (Optional)               #
##############################################
if args.align:
    print("Step 5: Running alignment...")
    if not args.starGenome or not args.geneAnnotation:
        print("Error: STAR genome directory and gene annotation file are required for alignment.")
        sys.exit(1)

    try:
        # Run STAR alignment
        print(f"Running STAR alignment with genome: {args.starGenome}")
        Splitseq_fun_lib.run_star_alignment_fun(args.numCores, args.starGenome, args.outputDir)
        print("STAR alignment completed.")

        # Run featureCounts
        print(f"Running featureCounts with annotation file: {args.geneAnnotation}")
        Splitseq_fun_lib.run_featureCounts_SAF_fun("GTF", args.numCores, args.geneAnnotation, args.outputDir)
        print("featureCounts completed.")

        # Run SAMtools
        print("Running SAMtools...")
        Splitseq_fun_lib.run_samtools_fun(args.outputDir)
        print("SAMtools processing completed.")

        # Run UMI-tools
        print("Running UMI-tools...")
        Splitseq_fun_lib.run_umi_tools_fun(args.outputDir)
        print("UMI-tools processing completed.")

    except Exception as e:
        print(f"Error during alignment or post-processing: {e}")
else:
    print("Step 5: Alignment step skipped.")

##############################################
# STEP 6: Final Cleanup                      #
##############################################
print("Step 6: Removing temporary files...")
for i in range(int(args.numCores)):
    Splitseq_fun_lib.remove_file_fun(f"split_fastq_F_{i}")
    Splitseq_fun_lib.remove_file_fun(f"split_fastq_R_{i}")

print("Pipeline completed successfully!")
