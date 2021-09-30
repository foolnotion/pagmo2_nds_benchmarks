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
/// Deductive sort
/**
 * Implementation of deductive sort. Complexity is \f$O M\kig N$\f (best) and \f$O MN^3$\f (worst).
 *
 * See: 
 *
 * McClymont and Keedwell 2012 - "Deductive Sort and Climbing Sort: New Methods for Non-Dominated Sorting"
 * Mishra S., Buzdalov M. 2020 - "If Unsure, Shuffle: Deductive Sort is Θ(MN3) but O(MN2) in Expectation over Input Permutations"
 *
 * @param points An std::vector containing the objectives of different individuals. Example
 * {{1,2,3},{-2,3,7},{-1,-2,-3},{0,0,0}}
 *
 * @return an std::tuple containing:
 *  - the non dominated fronts, an <tt>std::vector<std::vector<pop_size_t>></tt>
 * containing the non dominated fronts. Example {{1,2},{3},{0}}
 *  - the domination list, an <tt>std::vector<std::vector<pop_size_t>></tt>
 * containing the domination list, i.e. the indexes of all individuals
 * dominated by the individual at position \f$i\f$. Example {{},{},{0,3},{0}}
 *  - the domination count, an <tt>std::vector<pop_size_t></tt> containing the number of individuals
 * that dominate the individual at position \f$i\f$. Example {2, 0, 0, 1}
 *  - the non domination rank, an <tt>std::vector<pop_size_t></tt> containing the index of the non
 * dominated front to which the individual at position \f$i\f$ belongs. Example {2,0,0,1}
 *
 * @throws std::invalid_argument If the size of \p points is not at least 2
 *
 */
fnds_return_type deductive_sorting(const std::vector<vector_double> &points)
{
    auto N = points.size();
    // We make sure to have two points at least (one could also be allowed)
    if (N < 2u) {
        pagmo_throw(std::invalid_argument, "At least two points are needed for deductive_sorting: "
                                               + std::to_string(N) + " detected.");
    }

    std::vector<std::vector<pop_size_t>> non_dom_fronts;
    std::vector<std::vector<pop_size_t>> dom_list(N);
    std::vector<pop_size_t> dom_count(N);
    std::vector<pop_size_t> non_dom_rank(N);

    std::vector<bool> dominated(N, 0);
    std::vector<bool> sorted(N, 0);
    auto dominated_or_sorted = [&](size_t i) { return sorted[i] || dominated[i]; };

    decltype(N) n = 0;
    while (n < N) {
        std::vector<pop_size_t> front;

        for (decltype(N) i = 0; i < N; ++i) {
            if (!dominated_or_sorted(i)) {
                for (decltype(N) j = i + 1; j < N; ++j) {
                    if (!dominated_or_sorted(j)) {
                        dominated[i] = pareto_dominance(points[j], points[i]);
                        dominated[j] = !dominated[i] && pareto_dominance(points[i], points[j]); 

                        if (dominated[i]) {
                            dom_list[j].push_back(i);
                            ++dom_count[i];
                            break;
                        }
                    }
                }

                if (!dominated[i]) {
                    front.push_back(i);
                    non_dom_rank[i] = non_dom_fronts.size();
                    sorted[i] = 1;
                }
            }
        }
        std::fill(dominated.begin(), dominated.end(), 0ul);
        n += front.size();
        non_dom_fronts.push_back(front);
    }

    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}
} // namespace pagmo
