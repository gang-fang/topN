#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include "config.h"

namespace io {

/**
 * @brief Exception class for IO operations
 */
class IoException : public std::runtime_error {
public:
    explicit IoException(const std::string& message) 
        : std::runtime_error(message) {}
};

/**
 * @brief Class to handle FASTA format reading and writing
 */
class FastaReader {
public:
    explicit FastaReader(const std::string& filename);
    ~FastaReader();

    // Read all sequences from the file
    std::vector<Sequence> readAll();
    
    // Read next sequence (returns nullptr if end of file)
    std::unique_ptr<Sequence> readNext();
    
    // Reset to beginning of file
    void reset();

private:
    std::ifstream file;
    std::string currentLine;
    
    // Helper functions
    std::string readSequenceName();
    std::string readSequenceComment();
    std::string readSequenceData();
};

/**
 * @brief Class to write FASTA format
 */
class FastaWriter {
public:
    explicit FastaWriter(const std::string& filename);
    ~FastaWriter();

    // Write a single sequence
    void writeSequence(const Sequence& seq, int lineWidth = 60);
    
    // Write multiple sequences
    void writeSequences(const std::vector<Sequence>& sequences, int lineWidth = 60);
    
    // Write alignment results
    void writeResults(const AlignmentResults& results, bool includeAlignment = false);

private:
    std::ofstream file;
    void writeSequenceLine(const std::string& seq, int lineWidth);
};

// Utility functions
std::string trimWhitespace(const std::string& str);

} // namespace io

#endif // IO_UTILS_H 
