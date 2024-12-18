#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <sstream>

// Hamming distance function
int hammingDistance(const std::string &s1, const std::string &s2) {
    int distance = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) {
            ++distance;
        }
    }
    return distance;
}

// Filter barcode based on Hamming distance
std::string filterBarcode(const std::string &barcode, const std::unordered_set<std::string> &validBarcodes, int threshold) {
    for (const auto &valid : validBarcodes) {
        if (hammingDistance(barcode, valid) <= threshold) {
            return valid;
        }
    }
    return ""; // No match found
}

// Safe substring extraction
std::string safe_substr(const std::string &str, size_t start, size_t length) {
    if (start >= str.size()) return ""; // Start is out of bounds
    length = std::min(length, str.size() - start); // Adjust length to fit within bounds
    return str.substr(start, length);
}

// Learn barcode positions
void learnBarcodePositions(const std::string &fastqR, int &umi_start, int &umi_end,
                           int &bc1_start, int &bc1_end, int &bc2_start, int &bc2_end, int &bc3_start, int &bc3_end) {
    const std::string round1_static = "CATTCG";
    const std::string round2_static = "ATCCAC";
    const std::string round3_static = "GTGGCC";

    std::ifstream infile(fastqR);
    std::string line;
    int line_count = 0;

    std::vector<int> bc1_positions, bc2_positions, bc3_positions;

    while (std::getline(infile, line) && line_count < 4000) {
        if (line_count % 4 == 1) { // Sequence lines only
            if (line.find(round1_static) != std::string::npos)
                bc1_positions.push_back(line.find(round1_static));
            if (line.find(round2_static) != std::string::npos)
                bc2_positions.push_back(line.find(round2_static));
            if (line.find(round3_static) != std::string::npos)
                bc3_positions.push_back(line.find(round3_static));
        }
        line_count++;
    }

    // Default to valid bounds in case of missing sequences
    int seq_length = 101; // Adjust based on typical FASTQ sequence lengths
    umi_start = 0;
    umi_end = 10;
    bc1_start = 86;
    bc1_end = 94;
    bc2_start = 48;
    bc2_end = 56;
    bc3_start = 10;
    bc3_end = 18;

    // Validate positions if sequences were found
    if (!bc1_positions.empty() && !bc2_positions.empty() && !bc3_positions.empty()) {
        bc1_start = std::max(0, *std::max_element(bc1_positions.begin(), bc1_positions.end()) + 6);
        bc1_end = std::min(seq_length, bc1_start + 8);
        bc2_start = std::max(0, *std::max_element(bc2_positions.begin(), bc2_positions.end()) - 8);
        bc2_end = std::min(seq_length, bc2_start + 8);
        bc3_start = std::max(0, *std::max_element(bc3_positions.begin(), bc3_positions.end()) - 8);
        bc3_end = std::min(seq_length, bc3_start + 8);
        umi_start = std::max(0, bc3_start - 18);
        umi_end = std::min(seq_length, umi_start + 10);
    } else {
        std::cerr << "Warning: Static sequences not found in the provided fastqR file. Using default positions." << std::endl;
    }

    infile.close();
}

// Main demultiplexing function
extern "C" void demultiplex(const char *inputF, const char *inputR, const char *outputFile,
                            const char *round1File, const char *round2File, const char *round3File,
                            int errorThreshold, int numReadsBin) {
    std::unordered_set<std::string> round1Barcodes, round2Barcodes, round3Barcodes;

    // Load barcodes
    std::ifstream r1File(round1File), r2File(round2File), r3File(round3File);
    std::string barcode;
    while (r1File >> barcode) round1Barcodes.insert(barcode);
    while (r2File >> barcode) round2Barcodes.insert(barcode);
    while (r3File >> barcode) round3Barcodes.insert(barcode);

    r1File.close();
    r2File.close();
    r3File.close();

    // Learn barcode positions
    int umi_start, umi_end, bc1_start, bc1_end, bc2_start, bc2_end, bc3_start, bc3_end;
    learnBarcodePositions(inputR, umi_start, umi_end, bc1_start, bc1_end, bc2_start, bc2_end, bc3_start, bc3_end);

    // Process input FASTQs
    std::ifstream inF(inputF), inR(inputR);
    std::ofstream out(outputFile);
    std::string nameF, seqF, plusF, qualF;
    std::string nameR, seqR, plusR, qualR;

    int read_count = 0;
    while (std::getline(inF, nameF) && std::getline(inF, seqF) && std::getline(inF, plusF) && std::getline(inF, qualF) &&
           std::getline(inR, nameR) && std::getline(inR, seqR) && std::getline(inR, plusR) && std::getline(inR, qualR)) {
        if (read_count >= numReadsBin) break;

        // Extract UMI and barcodes safely
        std::string umi = safe_substr(seqR, umi_start, umi_end - umi_start);
        std::string barcode1 = safe_substr(seqR, bc1_start, bc1_end - bc1_start);
        std::string barcode2 = safe_substr(seqR, bc2_start, bc2_end - bc2_start);
        std::string barcode3 = safe_substr(seqR, bc3_start, bc3_end - bc3_start);

        // Validate barcodes
        std::string validBC1 = filterBarcode(barcode1, round1Barcodes, errorThreshold);
        std::string validBC2 = filterBarcode(barcode2, round2Barcodes, errorThreshold);
        std::string validBC3 = filterBarcode(barcode3, round3Barcodes, errorThreshold);

        // Add valid reads to the output
        if (!validBC1.empty() && !validBC2.empty() && !validBC3.empty()) {
            std::string newName = nameF + "_" + validBC1 + validBC2 + validBC3 + "_" + umi;
            out << newName << "\n" << seqF << "\n+\n" << qualF << "\n";
        }

        read_count++;
    }

    inF.close();
    inR.close();
    out.close();
}
