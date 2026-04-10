#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <vector>
#include <string>

struct Candidate {
    std::string targetName;
    double fastpScore;
    std::string alignedSegment;
};

using CandidateList = std::vector<Candidate>;
#endif
