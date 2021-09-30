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

#include <deque>

namespace pagmo
{
/// Hierarchical sorting
/**
 * An implementation of the efficient non-dominated sorting algorithm. Complexity is \f$O(M N \sqrt{N}$\f where \f$M\f$ is the number of objectives and \f$N\f$ is the number of individuals.
 * Worst-case complexity is \f$O(N^2)$\f.
 *
 * See: Bao et al. 2017 - "A novel non-dominated sorting algorithm for evolutionary multi-objective optimization"
 *      https://doi.org/10.1016/j.jocs.2017.09.015
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
 */
fnds_return_type hierarchical_sorting(const std::vector<vector_double> &points)
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

    std::deque<pop_size_t> q(N);
    std::iota(q.begin(), q.end(), 0ul);
    std::vector<pop_size_t> dominated;
    dominated.reserve(N);

    std::vector<std::vector<size_t>> fronts;
    while (!q.empty()) {
        std::vector<pop_size_t> front;

        while (!q.empty()) {
            auto q1 = q.front(); q.pop_front();
            front.push_back(q1);
            non_dom_rank[q1] = fronts.size();

            auto nonDominatedCount = 0ul;
            while (q.size() > nonDominatedCount) {
                auto qj = q.front(); q.pop_front();
                if (!pareto_dominance(points[q1], points[qj])) { 
                    q.push_back(qj);
                    ++nonDominatedCount;
                } else {
                    dominated.push_back(qj);
                }
            }
        }
        std::copy(dominated.begin(), dominated.end(), std::back_inserter(q));
        dominated.clear();
        fronts.push_back(front);

        std::stable_sort(q.begin(), q.end(),
            [&](auto i, auto j) {
                return std::lexicographical_compare(
                    points[i].begin(), points[i].end(),
                    points[j].begin(), points[j].end(),
                    std::ref(detail::less_than_f<vector_double::value_type>));
            });
    }
    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}
} // namespace pagmo
