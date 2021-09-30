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
/// Efficient non dominated sorting
/**
 * An implementation of the efficient non-dominated sorting algorithm. Complexity is \f$ O(MN \log N)\f$ (ENS-BS)
 * and \f$O(M N \sqrt{N}$\f (ENS-SS) where \f$M\f$ is the number of objectives and \f$N\f$ is the number of individuals.
 * Worst-case complexity is \f$O(N^2)$\f.
 *
 * See: Zhang et al. 2014 - "An Efficient Approach to Nondominated Sorting for Evolutionary Multiobjective Optimization"
 *      https://doi.org/10.1109/TEVC.2014.2308305
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
template<bool Binary>
inline fnds_return_type efficient_sorting(const std::vector<vector_double> &points)
{
    auto N = points.size();
    // We make sure to have two points at least (one could also be allowed)
    if (N < 2u) {
        pagmo_throw(std::invalid_argument, "At least two points are needed for fast_non_dominated_sorting: "
                                               + std::to_string(N) + " detected.");
    }
    // Initialize the return values
    std::vector<std::vector<pop_size_t>> non_dom_fronts;
    std::vector<std::vector<pop_size_t>> dom_list(N);
    std::vector<pop_size_t> dom_count(N);
    std::vector<pop_size_t> non_dom_rank(N);

    std::vector<pop_size_t> indices(N);
    std::iota(indices.begin(), indices.end(), 0ul);
    std::stable_sort(indices.begin(), indices.end(),
            [&](auto i, auto j) {
                return std::lexicographical_compare(
                    points[i].begin(), points[i].end(),
                    points[j].begin(), points[j].end(),
                    std::ref(detail::less_than_f<vector_double::value_type>));
            });

    // Check if individual i is dominated by any individual in the front f
    auto dominated = [&](const auto& front, auto i) {
        return std::any_of(front.rbegin(), front.rend(), [&](auto j) {
            return pareto_dominance(points[j], points[i]);
        });
    };

    for (auto i : indices) {
        decltype(non_dom_fronts)::iterator it;
        if constexpr (Binary) { // binary search
            it = std::partition_point(non_dom_fronts.begin(), non_dom_fronts.end(),
                        [&](const auto& front) { return dominated(front, i); });
        } else { // sequential search
            it = std::find_if(non_dom_fronts.begin(), non_dom_fronts.end(),
                        [&](const auto& front) { return !dominated(front, i); });
        }
        if (it == non_dom_fronts.end()) {
            non_dom_rank[i] = non_dom_fronts.size();
            non_dom_fronts.push_back({i});
        } else {
            non_dom_rank[i] = std::distance(non_dom_fronts.begin(), it);
            it->push_back(i);
        }
    }

    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}

fnds_return_type efficient_sorting_binary(const std::vector<vector_double> &points)
{
    return efficient_sorting<1>(points);
}

fnds_return_type efficient_sorting_sequential(const std::vector<vector_double> &points)
{
    return efficient_sorting<0>(points);
}

} // namespace pagmo
