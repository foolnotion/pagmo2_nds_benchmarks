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
// template<typename block_type = uint64_t, std::enable_if_t<std::is_unsigned_v<block_type>, bool> = true>
fnds_return_type rank_non_dominated_sorting(const std::vector<vector_double> &points)
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

    using block_type = uint64_t;
    using bitset = std::vector<block_type>; // a bitset is actually just a collection of 64-bit blocks
    size_t const block_bits = std::numeric_limits<block_type>::digits;
    size_t const num_blocks = (N / block_bits) + (N % block_bits != 0);
    block_type const block_max = std::numeric_limits<block_type>::max();

    std::vector<size_t> indices(N); // vector of indices that get sorted according to each objective
    std::iota(indices.begin(), indices.end(), 0ul);
    std::stable_sort(indices.begin(), indices.end(), [&](auto i, auto j) {
        return std::lexicographical_compare(points[i].begin(), points[i].end(), points[j].begin(), points[j].end(),
                                            std::ref(detail::less_than_f<vector_double::value_type>));
    });

    std::vector<block_type> b(num_blocks, block_max); // this bitset will help us compute the intersections
    b.back() >>= (block_bits * num_blocks - N);       // the bits beyond N in the last block must be set to zero

    std::vector<bitset> bs(N);                    // vector of bitsets (one for each individual)
    std::vector<std::pair<size_t, size_t>> br(N); // vector of ranges keeping track of the first/last non-zero blocks

    // compute initial bitsets
    for (auto i : indices) {
        b[i / block_bits] &= ~(block_type{1} << (i % block_bits)); // reset the bit corresponding to the current individual
        bs[i] = b;
        br[i] = {0, num_blocks - 1};
    }
    std::vector<bitset> rk(N);           // vector of sets keeping track of individuals whose rank was updated
    rk[0].resize(num_blocks, block_max); // initially all the individuals have the same rank

    // compute intersections
    for (decltype(M) k = 1; k < M; ++k) {
        std::stable_sort(indices.begin(), indices.end(),
                         [&](auto i, auto j) { return detail::less_than_f(points[i][k], points[j][k]); });
        std::fill(b.begin(), b.end(), block_max);
        b.back() >>= (block_bits * num_blocks - N);

        if (k < M - 1) {
            for (auto i : indices) {
                b[i / block_bits] &= ~(block_type{1} << (i % block_bits));
                auto [lo, hi] = br[i];
                if (lo > hi) {
                    continue;
                }

                auto p = bs[i].data();
                auto const q = b.data();

                // tighten the interval around non-zero blocks
                while (lo <= hi && !(p[lo] & q[lo]))
                    ++lo;
                while (lo <= hi && !(p[hi] & q[hi]))
                    --hi;

                // perform the blockwise intersection between sets
                for (size_t j = lo; j <= hi; ++j) {
                    p[j] &= q[j];
                }
                br[i] = {lo, hi};
            }
        } else {
            for (auto i : indices) {
                b[i / block_bits] &= ~(block_type{1} << (i % block_bits)); // reset bit for current individual
                auto [lo, hi] = br[i];
                if (lo > hi) {
                    continue;
                }

                auto p = bs[i].data();
                auto const q = b.data();

                // tighten the interval around non-zero blocks
                while (lo <= hi && !(p[lo] & q[lo]))
                    ++lo;
                while (lo <= hi && !(p[hi] & q[hi]))
                    --hi;

                auto r = rk[non_dom_rank[i]].data();
                auto rank_ = non_dom_rank[i] + 1;
                if (rk[rank_].empty()) {
                    rk[rank_].resize(num_blocks, 0);
                }

                // perform the blockwise intersection between sets
                for (size_t j = lo; j <= hi; ++j) {
                    auto o = block_bits * j;     // ordinal offset to get the individual index
                    auto v = p[j] & q[j] & r[j]; // obtain the dominance set
                    r[j] &= ~v;                  // remove dominance set individuals from current rank set

                    // iterate over the set bits in this block and update ranks
                    while (v) {
                        auto x = o + detail::count_trailing_zeros(v);
                        v &= (v - 1);
                        non_dom_rank[x] = rank_;
                        rk[rank_][x / block_bits] |= (block_type{1} << (x % block_bits));
                    }
                }
            }
        }
    }
    auto n_fronts = *std::max_element(non_dom_rank.begin(), non_dom_rank.end()) + 1;

    non_dom_fronts.resize(n_fronts);
    for (auto i : indices) {
        assert(i < non_dom_rank.size());
        assert(non_dom_rank[i] < non_dom_fronts.size());
        non_dom_fronts[non_dom_rank[i]].push_back(i);
    }
    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}
} // namespace pagmo
