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

#include <iostream>
#include <pagmo/utils/multi_objective.hpp>

namespace pagmo
{

namespace detail
{
using fnds_function_pointer_t = std::add_pointer_t<fnds_return_type(const std::vector<vector_double> &)>;

static const std::array<fnds_function_pointer_t, 7> algorithm_dispatch{
    &deductive_sorting,
    &efficient_sorting_binary,
    &efficient_sorting_sequential,
    &hierarchical_sorting,
    &merge_non_dominated_sorting,
    &fast_non_dominated_sorting,
    &rank_non_dominated_sorting
};
} // namespace detail

PAGMO_DLL_PUBLIC fnds_return_type non_dominated_sorting(const std::vector<vector_double> &points,
                                                        non_dominated_sorting_algorithm alg, bool dominate_on_equal)
{
    if (alg < non_dominated_sorting_algorithm::DS || alg > non_dominated_sorting_algorithm::RS) {
        pagmo_throw(std::invalid_argument, "unknown algorithm enum value: " + std::to_string(static_cast<int>(alg)));
    }

    std::vector<vector_double> work = points;
    std::stable_sort(work.begin(), work.end(), [](auto const& lhs, auto const& rhs) { 
        return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                                            std::ref(detail::less_than_f<vector_double::value_type>));
        });

    // handle duplicates
    // step 1: sort lexicographically
    std::vector<size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0ul);
    std::stable_sort(indices.begin(), indices.end(), [&](auto i, auto j) {
        return std::lexicographical_compare(points[i].begin(), points[i].end(), points[j].begin(), points[j].end(),
                                            std::ref(detail::less_than_f<vector_double::value_type>));
    });

    std::vector<std::pair<size_t, size_t>> duplicates;
    duplicates.reserve(points.size());

    // step 2: gather unique and duplicates
    // std::unique returns an iterator to the new logical end
    auto end = std::unique(indices.begin(), indices.end(), [&](auto i, auto j) {
        auto res = detail::equal_to_vf<vector_double::value_type>{}(points[i], points[j]);
        if (res) {
            duplicates.emplace_back(i, j);
        }
        return res;
    });

    std::vector<vector_double> unique;
    unique.reserve(points.size());

    std::vector<size_t> unique_idx;
    unique_idx.reserve(points.size());
    for (auto it = indices.begin(); it != end; ++it) {
        unique.push_back(points[*it]);
        unique_idx.push_back(*it);
    }

    // step 3: non-dominated sorting
    fnds_return_type retval = detail::algorithm_dispatch[static_cast<size_t>(alg)](unique);

    // these are the results with the individuals numbered from 0 to unique.size()
    // we need to map these results back to the original points vector
    auto const &[non_dom_fronts, dom_list, dom_count, non_dom_rank] = retval;

    std::vector<std::vector<pop_size_t>> fronts;
    std::vector<pop_size_t> rank(points.size());

    // translate ranks from indices in the unique vector to indices in the points vector
    size_t rank_max{0ul};
    for (decltype(unique.size()) i = 0; i < unique.size(); ++i) {
        size_t j = unique_idx[i];
        rank[j] = non_dom_rank[i];
        rank_max = std::max(rank_max, rank[j]);
    }

    // update ranks of duplicates
    for (auto [i, j] : duplicates) {
        rank[j] = rank[i] + dominate_on_equal;
        rank_max = std::max(rank_max, rank[j]);
    }

    // update fronts
    fronts.resize(rank_max + 1);
    for (size_t i = 0; i < rank.size(); ++i) {
        fronts[rank[i]].push_back(i);
    }

    // sort fronts (this ensures determinism across sorting algorithms which might otherwise
    // return the points in a front in a different order)
    for (auto &f : fronts) {
        std::sort(f.begin(), f.end());
    }

    // dom_list and dom_count are not complete but there's nothing we can do
    return std::make_tuple(std::move(fronts), std::move(dom_list), std::move(dom_count), std::move(rank));
}

} // namespace pagmo
