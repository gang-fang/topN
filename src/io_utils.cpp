#include "io_utils.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <sstream>
#include <iostream> // For std::cerr

namespace io {

// FastaReader implementation
FastaReader::FastaReader(const std::string& filename) {
    file.open(filename);
    if (!file) {
        throw IoException("Cannot open file: " + filename);
    }
}

FastaReader::~FastaReader() {
    if (file.is_open()) {
        file.close();
    }
}

std::vector<Sequence> FastaReader::readAll() {
    std::vector<Sequence> sequences;
    auto seq = readNext();
    while (seq) {
        sequences.push_back(*seq);
        seq = readNext();
    }
    return sequences;
}

std::unique_ptr<Sequence> FastaReader::readNext() {
    if (!file.is_open() || file.eof()) {
        return nullptr;
    }

    // Find next '>' character
    while (std::getline(file, currentLine)) {
        if (!currentLine.empty() && currentLine[0] == '>') {
            break;
        }
    }

    if (file.eof() && currentLine.empty()) {
        return nullptr;
    }

    // Parse the header line
    std::string name = readSequenceName();
    std::string comment = readSequenceComment();
    std::string sequence = readSequenceData();

    if (sequence.empty()) {
        throw IoException("Empty sequence for: " + name);
    }

    return std::make_unique<Sequence>(name, comment, sequence);
}

std::string FastaReader::readSequenceName() {
    if (currentLine.empty() || currentLine[0] != '>') {
        throw IoException("Invalid FASTA format: header line must start with '>'");
    }

    std::string header = currentLine.substr(1);  // Skip '>'
    std::string name;
    
    // Get first word as name
    std::istringstream iss(header);
    iss >> name;

    // Standard FASTA behavior for UniProt-style headers.
    size_t firstPipe = name.find('|');
    if (firstPipe != std::string::npos) {
        size_t secondPipe = name.find('|', firstPipe + 1);
        if (secondPipe != std::string::npos && secondPipe > firstPipe + 1) {
            name = name.substr(firstPipe + 1, secondPipe - firstPipe - 1);
        }
    }
    
    if (name.length() > MAX_SEQUENCE_NAME) {
        name = name.substr(0, MAX_SEQUENCE_NAME);
    }
    
    return name;
}

std::string FastaReader::readSequenceComment() {
    std::string header = currentLine.substr(1);  // Skip '>'
    
    // Find first space
    size_t spacePos = header.find_first_of(" \t");
    if (spacePos == std::string::npos) {
        return "";
    }
    
    std::string comment = header.substr(spacePos + 1);
    comment = trimWhitespace(comment);
    
    if (comment.length() > MAX_SEQUENCE_COMMENT) {
        comment = comment.substr(0, MAX_SEQUENCE_COMMENT);
    }
    
    return comment;
}

std::string FastaReader::readSequenceData() {
    std::string sequence;
    sequence.reserve(1024);  // Reserve some space

    // O(1) lookup table replacing the O(23) std::string::find per character
    static const auto makeAllowed = []() {
        std::array<bool, 256> t{};
        for (char c : {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X'})
            t[static_cast<unsigned char>(c)] = true;
        return t;
    };
    static const std::array<bool, 256> allowed = makeAllowed();

    while (std::getline(file, currentLine)) {
        if (currentLine.empty()) continue;
        if (currentLine[0] == '>') {
            file.seekg(-static_cast<long>(currentLine.length() + 1), std::ios::cur);
            break;
        }

        // Process each character: convert to uppercase and check if it's allowed
        for (char c : currentLine) {
            char upperC = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (allowed[static_cast<unsigned char>(upperC)]) {
                sequence += upperC;
            }
        }
    }

    return sequence;
}

void FastaReader::reset() {
    file.clear();
    file.seekg(0);
}

// FastaWriter implementation
FastaWriter::FastaWriter(const std::string& filename) {
    file.open(filename);
    if (!file) {
        throw IoException("Cannot open file for writing: " + filename);
    }
}

FastaWriter::~FastaWriter() {
    if (file.is_open()) {
        file.close();
    }
}

void FastaWriter::writeSequence(const Sequence& seq, int lineWidth) {
    file << '>' << seq.name;
    if (!seq.comment.empty()) {
        file << ' ' << seq.comment;
    }
    file << '\n';
    
    writeSequenceLine(seq.sequence, lineWidth);
}

void FastaWriter::writeSequences(const std::vector<Sequence>& sequences, int lineWidth) {
    for (const auto& seq : sequences) {
        writeSequence(seq, lineWidth);
    }
}

void FastaWriter::writeSequenceLine(const std::string& seq, int lineWidth) {
    for (size_t i = 0; i < seq.length(); i += lineWidth) {
        file << seq.substr(i, lineWidth) << '\n';
    }
}

void FastaWriter::writeResults(const AlignmentResults& results, bool includeAlignment) {
    for (const auto& result : results) {
        file << "Query: " << result.queryName << '\n'
             << "Target: " << result.targetName << '\n'
             << "Score: " << result.score << '\n'
             << "Identity: " << result.identity << '%' << '\n'
             << "Alignment length: " << result.alignmentLength << '\n';
        
        if (includeAlignment && !result.alignedQuery.empty()) {
            file << "\nAlignment:\n"
                 << result.alignedQuery << '\n'
                 << result.alignedTarget << '\n';
        }
        file << "\n";
    }
}

// Utility functions
std::string trimWhitespace(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

} // namespace io
