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

/// Non dominated front 2D (Kung's algorithm)
/**
 * Finds the non dominated front of a set of two dimensional objectives. Complexity is O(N logN) and is thus lower than
 * the
 * complexity of calling pagmo::fast_non_dominated_sorting
 *
 * See: Jensen, Mikkel T. "Reducing the run-time complexity of multiobjective EAs: The NSGA-II and other algorithms."
 * IEEE Transactions on Evolutionary Computation 7.5 (2003): 503-515.
 *
 * @param input_objs an <tt>std::vector</tt> containing the points (i.e. vector of objectives)
 *
 * @return A <tt>std::vector</tt> containing the indexes of the points in the non-dominated front
 *
 * @throws std::invalid_argument If the objective vectors are not all containing two-objectives
 */
std::vector<pop_size_t> non_dominated_front_2d(const std::vector<vector_double> &input_objs)
{
    // If the input is empty return an empty vector
    if (input_objs.size() == 0u) {
        return {};
    }
    // How many objectives? M, of course.
    auto M = input_objs[0].size();
    // We make sure all input_objs contain M objectives
    if (!std::all_of(input_objs.begin(), input_objs.end(),
                     [M](const vector_double &item) { return item.size() == M; })) {
        pagmo_throw(std::invalid_argument, "Input contains vector of objectives with heterogeneous dimensionalities");
    }
    // We make sure this function is only requested for two objectives.
    if (M != 2u) {
        pagmo_throw(std::invalid_argument, "The number of objectives detected is " + std::to_string(M)
                                               + ", while Kung's algorithm only works for two objectives.");
    }
    // Sanity checks are over. We may run Kung's algorithm.
    std::vector<pop_size_t> front;
    std::vector<pop_size_t> indexes(input_objs.size());
    std::iota(indexes.begin(), indexes.end(), pop_size_t(0u));
    // Sort in ascending order with respect to the first component
    std::sort(indexes.begin(), indexes.end(), [&input_objs](pop_size_t idx1, pop_size_t idx2) {
        if (detail::equal_to_f(input_objs[idx1][0], input_objs[idx2][0])) {
            return detail::less_than_f(input_objs[idx1][1], input_objs[idx2][1]);
        }
        return detail::less_than_f(input_objs[idx1][0], input_objs[idx2][0]);
    });
    for (auto i : indexes) {
        bool flag = false;
        for (auto j : front) {
            if (pareto_dominance(input_objs[j], input_objs[i])) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            front.push_back(i);
        }
    }
    return front;
}

/// Fast non dominated sorting
/**
 * An implementation of the fast non dominated sorting algorithm. Complexity is \f$ O(MN^2)\f$ where \f$M\f$ is the
 * number of objectives
 * and \f$N\f$ is the number of individuals.
 *
 * See: Deb, Kalyanmoy, et al. "A fast elitist non-dominated sorting genetic algorithm
 * for multi-objective optimization: NSGA-II." Parallel problem solving from nature PPSN VI. Springer Berlin Heidelberg,
 * 2000.
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
fnds_return_type fast_non_dominated_sorting(const std::vector<vector_double> &points)
{
    std::cout << "fast non-dominated sort\n";
    auto N = points.size();
    // We make sure to have two points at least (one could also be allowed)
    if (N < 2u) {
        pagmo_throw(std::invalid_argument, "At least two points are needed for fast_non_dominated_sorting: "
                                               + std::to_string(N) + " detected.");
    }
    // Initialize the return values
    std::vector<std::vector<pop_size_t>> non_dom_fronts(1u);
    std::vector<std::vector<pop_size_t>> dom_list(N);
    std::vector<pop_size_t> dom_count(N);
    std::vector<pop_size_t> non_dom_rank(N);

    // Start the fast non dominated sort algorithm
    for (decltype(N) i = 0u; i < N; ++i) {
        dom_list[i].clear();
        dom_count[i] = 0u;
        for (decltype(N) j = 0u; j < i; ++j) {
            if (pareto_dominance(points[i], points[j])) {
                dom_list[i].push_back(j);
                ++dom_count[j];
            } else if (pareto_dominance(points[j], points[i])) {
                dom_list[j].push_back(i);
                ++dom_count[i];
            }
        }
    }
    for (decltype(N) i = 0u; i < N; ++i) {
        if (dom_count[i] == 0u) {
            non_dom_rank[i] = 0u;
            non_dom_fronts[0].push_back(i);
        }
    }
    // we copy dom_count as we want to output its value at this point
    auto dom_count_copy(dom_count);
    auto current_front = non_dom_fronts[0];
    std::vector<std::vector<pop_size_t>>::size_type front_counter(0u);
    while (current_front.size() != 0u) {
        std::vector<pop_size_t> next_front;
        for (decltype(current_front.size()) p = 0u; p < current_front.size(); ++p) {
            for (decltype(dom_list[current_front[p]].size()) q = 0u; q < dom_list[current_front[p]].size(); ++q) {
                --dom_count_copy[dom_list[current_front[p]][q]];
                if (dom_count_copy[dom_list[current_front[p]][q]] == 0u) {
                    non_dom_rank[dom_list[current_front[p]][q]] = front_counter + 1u;
                    next_front.push_back(dom_list[current_front[p]][q]);
                }
            }
        }
        ++front_counter;
        current_front = next_front;
        if (current_front.size() != 0u) {
            non_dom_fronts.push_back(current_front);
        }
    }
    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}

} // namespace pagmo
