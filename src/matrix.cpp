#include "matrix.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>

bool loadSubstitutionMatrix(const std::string &filename,
                            std::vector<std::vector<int>> &matrix,
                            std::vector<char> &aminoAcids) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open substitution matrix file: " << filename << std::endl;
        return false;
    }
    return parseSubstitutionMatrix(infile, matrix, aminoAcids);
}

int getSubstitutionScore(char a, char b,
                         const std::vector<std::vector<int>> &matrix,
                         const std::vector<char> &aminoAcids) {
    int i = -1, j = -1;
    for (size_t idx = 0; idx < aminoAcids.size(); idx++) {
        if (std::toupper(a) == std::toupper(aminoAcids[idx]))
            i = idx;
        if (std::toupper(b) == std::toupper(aminoAcids[idx]))
            j = idx;
    }
    if (i == -1 || j == -1)
        return -1;
    return matrix[i][j];
} 
