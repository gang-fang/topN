#ifndef NW_ALIGN_H
#define NW_ALIGN_H

#include "config.h"  // for Sequence and ProgramConfig structs
#include <string>
#include <vector>
#include <utility>

/**
 * @brief Performs full-length global alignment with Hirschberg's divide-and-conquer traceback.
 *
 * @param query   The query sequence (as a std::string).
 * @param target  The target sequence (as a std::string).
 * @param config  Configuration settings (e.g., gap open and extend penalties).
 * @param freeEndGaps Reserved for API compatibility; currently ignored.
 * @return A pair of aligned sequences (first is the aligned query, second is the aligned target).
 */
std::pair<std::string, std::string> needlemanWunschAlignment(
    const std::string &query,
    const std::string &target,
    const ProgramConfig& config,
    bool freeEndGaps = false);

/**
 * @brief Computes the ranking similarity percentage with the distance-only DP used in the hot path.
 *
 * This does not build a traceback. It preprocesses both sequences and applies
 * the global distance formulation used for scoring candidate hits.
 *
 * @param querySeq  The query sequence.
 * @param targetSeq  The target sequence.
 * @param config    Configuration settings (e.g., gap open and extend penalties).
 * @param freeEndGaps Reserved for API compatibility; currently ignored.
 * @return The similarity percentage (0.0-100.0).
 */
double computeNWSimilarityPercentage(
    const Sequence &querySeq,
    const Sequence &targetSeq,
    const ProgramConfig& config,
    bool freeEndGaps = false);

#endif  // NW_ALIGN_H
