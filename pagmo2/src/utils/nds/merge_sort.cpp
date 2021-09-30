/* Copyright 2017-2021 PaGMO development team

This file is part of the PaGMO library.

The PaGMO library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The PaGMO library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the PaGMO library.  If not,
see https://www.gnu.org/licenses/. */

#include <pagmo/utils/multi_objective.hpp>

namespace pagmo
{
namespace detail
{
constexpr int INSERTIONSORT = 7;

class BitsetManager
{
    using word_t = uint64_t;

    static constexpr int FIRST_WORD_RANGE = 0;
    static constexpr int LAST_WORD_RANGE = 1;
    static constexpr int N_BIT_ADDR = 6; // 2^6 = 64
    static constexpr word_t WORD_MASK = ~word_t{0};
    static constexpr int WORD_SIZE = std::numeric_limits<word_t>::digits;

    std::vector<std::vector<word_t>> bitsets;
    std::vector<std::array<int, 2>> bsRanges;
    std::vector<pop_size_t> wordRanking; // Ranking of each bitset word. A bitset word contains 64 solutions.
    std::vector<pop_size_t> ranking, ranking0;
    int maxRank = 0;
    std::vector<long> incrementalBitset;
    int incBsFstWord, incBsLstWord;

public:
    std::vector<pop_size_t> const &getRanking() const
    {
        return ranking0;
    }

    bool updateSolutionDominance(int solutionId)
    {
        int fw = bsRanges[solutionId][FIRST_WORD_RANGE];
        int lw = bsRanges[solutionId][LAST_WORD_RANGE];
        if (lw > incBsLstWord) {
            lw = incBsLstWord;
        }
        if (fw < incBsFstWord) {
            fw = incBsFstWord;
        }

        while (fw <= lw && 0 == (bitsets[solutionId][fw] & incrementalBitset[fw])) {
            fw++;
        }
        while (fw <= lw && 0 == (bitsets[solutionId][lw] & incrementalBitset[lw])) {
            lw--;
        }
        bsRanges[solutionId][FIRST_WORD_RANGE] = fw;
        bsRanges[solutionId][LAST_WORD_RANGE] = lw;

        if (fw > lw) {
            return false;
        }
        for (; fw <= lw; fw++) {
            bitsets[solutionId][fw] &= incrementalBitset[fw];
        }
        return true;
    }

    void computeSolutionRanking(int solutionId, int initSolId)
    {
        int fw = bsRanges[solutionId][FIRST_WORD_RANGE];
        int lw = bsRanges[solutionId][LAST_WORD_RANGE];
        if (lw > incBsLstWord) {
            lw = incBsLstWord;
        }
        if (fw < incBsFstWord) {
            fw = incBsFstWord;
        }
        if (fw > lw) {
            return;
        }
        word_t word;
        int i = 0, rank = 0, offset;

        for (; fw <= lw; fw++) {
            word = bitsets[solutionId][fw] & incrementalBitset[fw];
            if (word != 0) {
                i = (int)detail::count_trailing_zeros(static_cast<word_t>(word));
                offset = fw * WORD_SIZE;
                do {
                    if (ranking[offset + i] >= rank) {
                        rank = ranking[offset + i] + 1;
                    }
                    i++;
                    word_t w = static_cast<word_t>(word) >> i;
                    i += w ? (int)detail::count_trailing_zeros(w) : (int)WORD_SIZE;
                } while (i < WORD_SIZE && rank <= wordRanking[fw]);
                if (rank > maxRank) {
                    maxRank = rank;
                    break;
                }
            }
        }
        ranking[solutionId] = rank;
        ranking0[initSolId] = rank;
        i = solutionId >> N_BIT_ADDR;
        if (rank > wordRanking[i]) {
            wordRanking[i] = rank;
        }
    }

    void updateIncrementalBitset(int solutionId)
    {
        int wordIndex = solutionId >> N_BIT_ADDR;
        // int shiftDistance = solutionId & 0x3f;
        int shiftDistance = solutionId;
        incrementalBitset[wordIndex] |= (word_t{1} << shiftDistance);
        if (incBsLstWord < wordIndex) incBsLstWord = wordIndex;
        if (incBsFstWord > wordIndex) incBsFstWord = wordIndex;
    }

    bool initializeSolutionBitset(int solutionId)
    {
        // int const shiftDistance = solutionId & 0x3f;
        int const shiftDistance = solutionId;
        int wordIndex = solutionId >> N_BIT_ADDR;
        if (wordIndex < incBsFstWord || 0 == solutionId) {
            bsRanges[solutionId][FIRST_WORD_RANGE] = std::numeric_limits<int>::max();
            return false;
        } else if (wordIndex == incBsFstWord) { // only 1 word in common
            bitsets[solutionId] = std::vector<word_t>(wordIndex + 1);
            long intersection = incrementalBitset[incBsFstWord] & ~(WORD_MASK << shiftDistance);
            if (intersection != 0) {
                bsRanges[solutionId][FIRST_WORD_RANGE] = wordIndex;
                bsRanges[solutionId][LAST_WORD_RANGE] = wordIndex;
                bitsets[solutionId][wordIndex] = intersection;
            }
            return intersection != 0;
        }
        // more than one word in common
        int lw = incBsLstWord < wordIndex ? incBsLstWord : wordIndex;
        bsRanges[solutionId][FIRST_WORD_RANGE] = incBsFstWord;
        bsRanges[solutionId][LAST_WORD_RANGE] = lw;
        bitsets[solutionId] = std::vector<word_t>(lw + 1);
        std::copy_n(incrementalBitset.begin() + incBsFstWord, lw - incBsFstWord + 1,
                    bitsets[solutionId].begin() + incBsFstWord);
        if (incBsLstWord >= wordIndex) { // update (compute intersection) the last word
            bitsets[solutionId][lw] = incrementalBitset[lw] & ~(WORD_MASK << shiftDistance);
            if (bitsets[solutionId][lw] == 0) {
                bsRanges[solutionId][LAST_WORD_RANGE]--;
            }
        }
        return true;
    }

    void clearIncrementalBitset()
    {
        std::fill(incrementalBitset.begin(), incrementalBitset.end(), 0ul);
        incBsLstWord = 0;
        incBsFstWord = std::numeric_limits<int>::max();
        maxRank = 0;
    }

    BitsetManager() = default;

    // constructor
    BitsetManager(size_t nSolutions)
    {
        int n = (int)nSolutions - 1;
        int wordIndex = static_cast<int>(static_cast<size_t>(n) >> N_BIT_ADDR);
        ranking.resize(nSolutions, 0);
        ranking0.resize(nSolutions, 0);
        wordRanking.resize(nSolutions, 0);
        bitsets.resize(nSolutions);
        bsRanges.resize(nSolutions);
        incrementalBitset.resize(wordIndex + 1);
        incBsLstWord = 0;
        incBsFstWord = std::numeric_limits<int>::max();
        maxRank = 0;
    }
};

inline int CompareLex(vector_double const &s1, vector_double const &s2, size_t fromObj, size_t toObj)
{
    for (; fromObj < toObj; fromObj++) {
        if (s1[fromObj] < s2[fromObj]) return -1;
        if (s1[fromObj] > s2[fromObj]) return 1;
    }
    return 0;
}

inline bool MergeSort(std::vector<vector_double> &src, std::vector<vector_double> &dest, int low, int high, int obj,
                      int toObj)
{
    int i, j, s;
    int destLow = low;
    int length = high - low;

    if (length < INSERTIONSORT) {
        bool alreadySorted{true};
        for (i = low; i < high; i++) {
            for (j = i; j > low && CompareLex(dest[j - 1], dest[j], obj, toObj) > 0; j--) {
                alreadySorted = false;
                dest[j].swap(dest[j - 1]);
            }
        }
        return alreadySorted; // if temp==null, src is already sorted
    }
    int mid = (low + high) / 2;
    bool isSorted = MergeSort(dest, src, low, mid, obj, toObj) & MergeSort(dest, src, mid, high, obj, toObj);

    // If list is already sorted, just copy from src to dest.
    if (src[mid - 1][obj] <= src[mid][obj]) {
        std::copy_n(src.begin() + low, length, dest.begin() + destLow);
        return isSorted;
    }

    for (s = low, i = low, j = mid; s < high; s++) {
        if (j >= high) {
            dest[s] = src[i++];
        } else if (i < mid && CompareLex(src[i], src[j], obj, toObj) <= 0) {
            dest[s] = src[i++];
        } else {
            dest[s] = src[j++];
        }
    }
    return false;
}
} // namespace detail

fnds_return_type merge_non_dominated_sorting(const std::vector<vector_double> &points)
{
    auto N = points.size();
    auto M = points.front().size();
    // We make sure to have two points at least (one could also be allowed)
    if (N < 2u) {
        pagmo_throw(std::invalid_argument,
                    "At least two points are needed for non-dominated sorting: " + std::to_string(N) + " detected.");
    }

    // Initialize the return values
    std::vector<std::vector<pop_size_t>> non_dom_fronts;
    std::vector<std::vector<pop_size_t>> dom_list(N);
    std::vector<pop_size_t> dom_count(N);
    std::vector<pop_size_t> non_dom_rank(N, 0);

    detail::BitsetManager bsm(N);
    std::vector<vector_double> population(N);
    std::vector<vector_double> work(N);

    const pop_size_t SOL_ID = M;
    const pop_size_t SORT_INDEX = SOL_ID + 1;

    // initialize auxiliary data structures
    for (decltype(N) i = 0; i < N; ++i) {
        population[i].resize(SORT_INDEX + 1);
        std::copy_n(points[i].begin(), M, population[i].begin());
        population[i][SOL_ID] = static_cast<double>(i);
        population[i][SORT_INDEX] = static_cast<double>(i); // points are already sorted
    }

    decltype(N) solutionId;
    bool dominance{false};
    work = population;
    detail::MergeSort(population, work, 0, N, 1, 2);
    population = work;
    for (decltype(N) p = 0; p < N; p++) {
        solutionId = (int)population[p][SORT_INDEX];
        dominance |= bsm.initializeSolutionBitset(solutionId);
        bsm.updateIncrementalBitset(solutionId);
        if (2 == M) {
            int initSolId = (int)population[p][SOL_ID];
            bsm.computeSolutionRanking(solutionId, initSolId);
        }
    }

    if (M > 2) {
        dominance = false;
        decltype(M) lastObjective = M - 1;
        work = population;
        for (int obj = 2; obj < M; obj++) {
            if (detail::MergeSort(population, work, 0, N, obj,
                                  obj + 1)) { // Population has the same order as in previous objective
                if (obj == lastObjective) {
                    for (decltype(N) p = 0; p < N; p++)
                        bsm.computeSolutionRanking((int)population[p][SORT_INDEX], (int)population[p][SOL_ID]);
                }
                continue;
            }
            population = work;
            bsm.clearIncrementalBitset();
            dominance = false;
            for (decltype(N) p = 0; p < N; p++) {
                auto initSolId = population[p][SOL_ID];
                solutionId = population[p][SORT_INDEX];
                if (obj < lastObjective) {
                    dominance |= bsm.updateSolutionDominance(solutionId);
                } else {
                    bsm.computeSolutionRanking(solutionId, initSolId);
                }
                bsm.updateIncrementalBitset(solutionId);
            }
        }
    }

    non_dom_rank = bsm.getRanking();

    non_dom_fronts.resize(*std::max_element(non_dom_rank.begin(), non_dom_rank.end()) + 1);
    for (decltype(N) i = 0; i < N; i++) {
        non_dom_fronts[non_dom_rank[i]].push_back(i);
    }
    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
            std::move(non_dom_rank));
}

} // namespace pagmo
