import ctypes
import os

# Load the C++ shared library
lib_path = "./barcode_demux.so"  # Change path for Windows if needed
barcode_demux = ctypes.CDLL(lib_path)

# Define the argument types and return type for the C++ function
barcode_demux.demultiplex.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
barcode_demux.demultiplex.restype = None

def run_demultiplex(inputF, inputR, output, barcode_file, error_threshold):
    """
    Wrapper for the C++ demultiplex function.
    """
    barcode_demux.demultiplex(inputF.encode('utf-8'), inputR.encode('utf-8'),
                              output.encode('utf-8'), barcode_file.encode('utf-8'),
                              error_threshold)

# Example usage
if __name__ == "__main__":
    input_fastqF = "split_fastq_F_0.fastq"
    input_fastqR = "split_fastq_R_0.fastq"
    output_file = "output_demux.fastq"
    barcode_file = "Round1_barcodes_new5.txt"
    error_threshold = 1

    run_demultiplex(input_fastqF, input_fastqR, output_file, barcode_file, error_threshold)
    print("Demultiplexing completed.")
