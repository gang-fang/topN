#ifndef FASTP_H
#define FASTP_H

#include "candidate.h"  // Access Candidate and CandidateList
#include "config.h"     // For Sequence and ProgramConfig structs
#include <string>
#include <vector>
#include <unordered_map>

/**
 * @brief Computes a fast prefilter score between query and target protein sequences.
 *
 * For each k-mer in the query, the function checks if an exact k-mer match appears in the target.
 * If so, a bonus is added to the score. Uses HSP scoring, normalized by sequence length.
 *
 * @param query The query protein sequence.
 * @param target The target protein sequence.
 * @param config Configuration settings (e.g., k-mer size, lower threshold).
 * @param aminoAcids The amino acid symbols corresponding to the matrix.
 * @param alignedSegment Optional pointer to store the aligned segment.
 * @return A float value representing the computed similarity.
 */
float computeFastpScore(const Sequence& query,
                        const Sequence& target,
                        const ProgramConfig& config,
                        const std::vector<char>& aminoAcids,
                        std::string* alignedSegment = nullptr);

/**
 * @brief Runs the FASTP prefilter on a set of query and database sequences.
 *
 * For every query, it computes a fastp score against each database protein and then
 * selects the top candidates based on the configuration.
 *
 * @param queries Vector of query sequences.
 * @param db Vector of database (target) sequences.
 * @param config Configuration settings (e.g., topK, k-mer size, lower threshold).
 * @param aminoAcids The corresponding amino acid symbols.
 * @return A vector of CandidateList (one list per query).
 */
std::vector<CandidateList> runFastpPrefilter(const std::vector<Sequence>& queries,
                                             const std::vector<Sequence>& db,
                                             const ProgramConfig& config,
                                             const std::vector<char>& aminoAcids);

/**
 * @brief Checks if the length ratio between two sequences is acceptable.
 * 
 * @param query The query sequence.
 * @param target The target sequence.
 * @param lengthRatio The maximum acceptable length ratio between sequences.
 * @return true if the ratio is less than or equal to the specified lengthRatio.
 */
bool isLengthRatioAcceptable(const Sequence& query, const Sequence& target, float lengthRatio);

#endif  // FASTP_H
