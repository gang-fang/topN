#include "nw_align.h"
#include "speedup.h"
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cctype>
#include <climits>
#include <cstdint>

namespace {
    // -----------------------------------------------------------------------
    // Hirschberg divide-and-conquer Needleman-Wunsch traceback used only by
    // needlemanWunschAlignment(); the ranking hot path uses the distance-only
    // DP further below.
    // -----------------------------------------------------------------------

    class BufferPool {
    public:
        void ensureSize(size_t size) {
            if (size > bufferSize) {
                ccArray = std::make_unique<int[]>(size);
                ddArray = std::make_unique<int[]>(size);
                rrArray = std::make_unique<int[]>(size);
                ssArray = std::make_unique<int[]>(size);
                bufferSize = size;
            }
        }
        int* cc() { return ccArray.get(); }
        int* dd() { return ddArray.get(); }
        int* rr() { return rrArray.get(); }
        int* ss() { return ssArray.get(); }
    private:
        std::unique_ptr<int[]> ccArray, ddArray, rrArray, ssArray;
        size_t bufferSize = 0;
    };

    static thread_local BufferPool pool;

    enum class Op { Replace = 0, Delete, Insert };

    class Script {
    public:
        void append(Op op, int count = 1) {
            if (op == Op::Replace) {
                ops.push_back(0);
            } else if (op == Op::Delete) {
                if (!ops.empty() && ops.back() < 0) ops.back() -= count;
                else ops.push_back(-count);
            } else {
                if (!ops.empty() && ops.back() > 0) ops.back() += count;
                else ops.push_back(count);
            }
        }
        const std::vector<int>& operations() const { return ops; }
    private:
        std::vector<int> ops;
    };

    class Aligner {
    public:
        Aligner(const std::string& seqA, const std::string& seqB, const ProgramConfig& config)
            : A(seqA), B(seqB),
              g(config.gapOpenPenalty - config.gapExtendPenalty),
              h(config.gapExtendPenalty),
              m(config.gapOpenPenalty + config.gapExtendPenalty)
        {}

        std::pair<std::string, std::string> align() {
            script = std::make_unique<Script>();
            score = alignRecursive(0, A.length(), 0, B.length(), -g, -g);
            return buildAlignment();
        }

        int getScore() const { return score; }

    private:
        const std::string& A, &B;
        int g, h, m;
        int score = 0;
        std::unique_ptr<Script> script;

        int gapCost(int k) const { return k <= 0 ? 0 : g + h * k; }

        int alignRecursive(int startA, int endA, int startB, int endB, int tb, int te) {
            int lenA = endA - startA, lenB = endB - startB;
            if (lenB <= 0) {
                if (lenA > 0) script->append(Op::Delete, lenA);
                return -gapCost(lenA);
            }
            if (lenA <= 1) {
                if (lenA <= 0) { script->append(Op::Insert, lenB); return -gapCost(lenB); }
                return alignSingleA(startA, startB, endB, tb, te);
            }
            int midi = startA + lenA / 2;
            auto [midj, midc, type] = findMidpoint(startA, midi, endA, startB, endB, tb, te);
            if (type == 1) {
                return alignRecursive(startA, midi, startB, midj, tb, -g)
                     + alignRecursive(midi, endA, midj, endB, -g, te);
            } else {
                alignRecursive(startA, midi - 1, startB, midj, tb, 0);
                script->append(Op::Delete, 2);
                alignRecursive(midi + 1, endA, midj, endB, 0, te);
                return midc;
            }
        }

        int alignSingleA(int posA, int startB, int endB, int tb, int te) {
            int lenB = endB - startB;
            int midc = (tb < te ? te : tb) - h - gapCost(lenB);
            int midj = 0;
            for (int j = 0; j < lenB; ++j) {
                int c = -gapCost(j) + substitutionMatrix[A[posA] - 'A'][B[startB + j] - 'A'] - gapCost(lenB - j - 1);
                if (c > midc) { midc = c; midj = j + 1; }
            }
            if (midj == 0) {
                script->append(Op::Insert, lenB);
                script->append(Op::Delete, 1);
            } else {
                if (midj > 1) script->append(Op::Insert, midj - 1);
                script->append(Op::Replace);
                if (midj < lenB) script->append(Op::Insert, lenB - midj);
            }
            return midc;
        }

        std::tuple<int, int, int> findMidpoint(int startA, int midi, int endA,
                                                int startB, int endB, int tb, int te) {
            int lenB = endB - startB;
            pool.ensureSize(lenB + 1);
            int* cc = pool.cc(); int* dd = pool.dd();
            int* rr = pool.rr(); int* ss = pool.ss();

            cc[0] = 0;
            int t = -g;
            for (int j = 1; j <= lenB; ++j) { cc[j] = t = t - h; dd[j] = t - g; }
            t = tb;
            for (int i = startA + 1; i <= midi; ++i) {
                int s = cc[0]; cc[0] = t = t - h; int e = t - g;
                for (int j = 1; j <= lenB; ++j) {
                    int c, d;
                    e = e - h; if ((c = cc[j] - m) > e) e = c;
                    d = dd[j] - h; if ((c = cc[j] - m) > d) d = c;
                    c = s + substitutionMatrix[A[i-1]-'A'][B[startB+j-1]-'A'];
                    if (e > c) c = e;
                    if (d > c) c = d;
                    s = cc[j]; cc[j] = c; dd[j] = d;
                }
            }
            dd[0] = cc[0];

            rr[lenB] = 0; t = -g;
            for (int j = lenB - 1; j >= 0; --j) { rr[j] = t = t - h; ss[j] = t - g; }
            t = te;
            for (int i = endA; i > midi; --i) {
                int s = rr[lenB]; rr[lenB] = t = t - h; int e = t - g;
                for (int j = lenB - 1; j >= 0; --j) {
                    int c, d;
                    e = e - h; if ((c = rr[j] - m) > e) e = c;
                    d = ss[j] - h; if ((c = rr[j] - m) > d) d = c;
                    c = s + substitutionMatrix[A[i-1]-'A'][B[startB+j]-'A'];
                    if (e > c) c = e;
                    if (d > c) c = d;
                    s = rr[j]; rr[j] = c; ss[j] = d;
                }
            }
            ss[lenB] = rr[lenB];

            int midc = cc[0] + rr[0], midj = 0, type = 1;
            for (int j = 0; j <= lenB; ++j) {
                int c = cc[j] + rr[j];
                if (c >= midc && (c > midc || (cc[j] != dd[j] && rr[j] == ss[j])))
                    { midc = c; midj = j; }
            }
            for (int j = 0; j <= lenB; ++j) {
                int c = dd[j] + ss[j] + g;
                if (c > midc) { midc = c; midj = j; type = 2; }
            }
            if (startA == 0 && static_cast<size_t>(endA) == A.length() &&
                startB == 0 && static_cast<size_t>(endB) == B.length())
                score = midc;
            return {startB + midj, midc, type};
        }

        std::pair<std::string, std::string> buildAlignment() {
            std::string alignedA, alignedB;
            int i = 0, j = 0;
            for (int op : script->operations()) {
                if (op == 0) { alignedA += A[i++]; alignedB += B[j++]; }
                else if (op > 0) { for (int k=0;k<op;++k){alignedA+='-';alignedB+=B[j++];} }
                else { for (int k=0;k<-op;++k){alignedA+=A[i++];alignedB+='-';} }
            }
            return {alignedA, alignedB};
        }
    };

    static thread_local std::vector<int> distH;
    static thread_local std::vector<int> distF;
    static thread_local std::vector<int> distE;

    constexpr int DIST_INF = INT_MAX / 4;

    double computeGlobalSimilarity(const std::string& A, const std::string& B,
                                         int open, int ext) {
        const int lenA = static_cast<int>(A.size());
        const int lenB = static_cast<int>(B.size());
        const int cols = lenB + 1;

        if (static_cast<int>(distH.size()) < cols) {
            distH.resize(cols);
            distF.resize(cols);
            distE.resize(cols);
        }

        int* H = distH.data();
        int* F = distF.data();
        int* E = distE.data();

        H[0] = 0;
        E[0] = DIST_INF;
        F[0] = DIST_INF;
        for (int j = 1; j <= lenB; ++j) {
            H[j] = open + (j - 1) * ext;
            E[j] = H[j];
            F[j] = DIST_INF;
        }

        for (int i = 1; i <= lenA; ++i) {
            int prevDiag = H[0];
            H[0] = open + (i - 1) * ext;
            F[0] = H[0];
            E[0] = DIST_INF;

            const unsigned char ai = static_cast<unsigned char>(A[i - 1]);
            int e = DIST_INF;
            for (int j = 1; j <= lenB; ++j) {
                const int up = H[j];
                e = std::min(H[j - 1] + open, e + ext);
                F[j] = std::min(up + open, F[j] + ext);

                const unsigned char bj = static_cast<unsigned char>(B[j - 1]);
                const int sub = substitutionWeightMatrix[ai - 'A'][bj - 'A'];
                H[j] = std::min({prevDiag + sub, e, F[j]});
                prevDiag = up;
            }
        }

        const int norm = std::max(lenA, lenB);
        return norm > 0 ? 100.0 - static_cast<double>(H[lenB]) / static_cast<double>(norm) : 0.0;
    }

} // anonymous namespace

// Normalize sequences before alignment/scoring by uppercasing and dropping
// non-alphabetic characters.
std::string preprocess(const std::string& seq) {
    std::string result;
    result.reserve(seq.length());
    for (char c : seq)
        if (std::isalpha(c)) result += std::toupper(c);
    return result;
}

// Full traceback-producing alignment for the public API. This is separate from
// the scoring-only DP used during ranking.
std::pair<std::string, std::string> needlemanWunschAlignment(
    const std::string& query, const std::string& target,
    const ProgramConfig& config, [[maybe_unused]] bool freeEndGaps) {
    std::string seqA = preprocess(query);
    std::string seqB = preprocess(target);
    if (seqA.empty() || seqB.empty())
        throw std::invalid_argument("Empty sequence after preprocessing");
    Aligner aligner(seqA, seqB, config);
    return aligner.align();
}

// Similarity percentage for ranking, using the global distance
// formulation with normalized substitution weights in [0,100], affine gaps
// k0 + (l-1) * k1, and final score 100 - distance / max(lenA, lenB).
double computeNWSimilarityPercentage(
    const Sequence& querySeq, const Sequence& targetSeq,
    const ProgramConfig& config, [[maybe_unused]] bool freeEndGaps) {
    std::string seqA = preprocess(querySeq.sequence);
    std::string seqB = preprocess(targetSeq.sequence);
    if (seqA.empty() || seqB.empty()) return 0.0;
    return computeGlobalSimilarity(seqA, seqB,
                                         config.gapOpenPenalty,
                                         config.gapExtendPenalty);
}
