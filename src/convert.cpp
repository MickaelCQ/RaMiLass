//
// Created by raphael on 11/13/25.
//

#include "convert.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <algorithm> // For std::remove_if
#include <cctype>    // For ::isspace

void Convert::process_fasta_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    // Clear any previous data
    bit_vector.clear();
    read_end_positions.clear();

    size_t total_read_size = 0;
    size_t total_read_number = 0;
    std::string line;

    // --- FIRST PASS: Calculate total size and read count ---
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }

        if (line[0] == '>') {
            total_read_number++; // Found a new read header
        } else {
            // This is a sequence line.
            // Remove any potential whitespace from the line before counting.
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            total_read_size += line.length();
        }
    }

    // --- MEMORY RESERVATION ---
    // Reserve space for all bits (2 bits per nucleotide)
    bit_vector.reserve(total_read_size * 2);
    // Reserve space for all read end positions
    read_end_positions.reserve(total_read_number);

    // --- SECOND PASS: Convert and store data ---
    // Reset file stream to the beginning
    file.clear();
    file.seekg(0, std::ios::beg);

    std::stringstream current_sequence_stream;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }

        if (line[0] == '>') {
            // This is a header line.
            // Process the sequence we've built up so far.
            std::string sequence = current_sequence_stream.str();
            if (!sequence.empty()) {
                convert_and_store_sequence(sequence);
            }

            // Reset the stream for the new sequence
            current_sequence_stream.str("");
            current_sequence_stream.clear();
        } else {
            // This is a sequence line. Append it to our current sequence stream.
            // We strip any potential whitespace from the line.
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            current_sequence_stream << line;
        }
    }

    // After the loop, process the very last sequence in the file
    std::string last_sequence = current_sequence_stream.str();
    if (!last_sequence.empty()) {
        convert_and_store_sequence(last_sequence);
    }

    file.close();
}

/**
 * @brief (Private Helper) Converts a single DNA sequence string and appends it to the bit vector.
 *
 * This helper function iterates through the sequence, converts each nucleotide
 * to two bits, and appends them to the main `bit_vector`. It then records
 * the new total size of `bit_vector` in `read_end_positions`.
 *
 * This function assumes `bit_vector` has already been reserved.
 */
void Convert::convert_and_store_sequence(const std::string& sequence) {
    if (sequence.empty()) {
        return;
    }

    // NOTE: We no longer call bit_vector.reserve() here,
    // as it was done in process_fasta_file().

    for (char nucleotide : sequence) {
        // Handle both upper and lower case
        switch (toupper(nucleotide)) {
            case 'A':
                bit_vector.push_back(false); // 0
                bit_vector.push_back(false); // 0
                break;
            case 'C':
                bit_vector.push_back(true);  // 1
                bit_vector.push_back(false); // 0
                break;
            case 'G':
                bit_vector.push_back(false); // 0
                bit_vector.push_back(true);  // 1
                break;
            case 'T':
                bit_vector.push_back(true);  // 1
                bit_vector.push_back(true);  // 1
                break;
            default:
                // Handle N or other non-standard characters.
                // We will treat them as 'A' (00) for this example,
                // but you can skip them or throw an error.
                // Skipping them would misalign the bit vector and read positions.
                // Let's treat them as 'A' (00) to maintain sequence length.
                if (nucleotide != '\n' && nucleotide != '\r') {
                     std::cerr << "Warning: Unhandled character '" << nucleotide << "'. Treating as 'A'." << std::endl;
                }
                bit_vector.push_back(false); // 0
                bit_vector.push_back(false); // 0
                break;
        }
    }

    // Store the *new* cumulative size of the bit vector
    // This marks the end position of the read we just added.
    read_end_positions.push_back(bit_vector.size());
}

const std::vector<bool>& Convert::get_bit_vector() const {
    return bit_vector;
}

const std::vector<size_t>& Convert::get_read_end_positions() const {
    return read_end_positions;
}