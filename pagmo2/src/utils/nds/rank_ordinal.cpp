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

#include <Eigen/Core>

namespace pagmo
{

// template<typename block_type = uint64_t, std::enable_if_t<std::is_unsigned_v<block_type>, bool> = true>
fnds_return_type rank_ordinal_non_dominated_sorting(const std::vector<vector_double> &points)
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

    using Mat = Eigen::Array<size_t, -1, -1>;
    using Vec = Eigen::Array<size_t, 1, -1>;

    // 1) sort indices according to the stable sorting rules
    // already sorted!
    //std::stable_sort(indices.begin(), indices.end(), [&](auto const i, auto const j) {
    //    return std::lexicographical_compare(points[i].begin(), points[i].end(), points[j].begin(), points[j].end(),
    //                                        std::ref(detail::less_than_f<vector_double::value_type>));
    //});
    Mat p(N, M); // permutation matrix
    Mat r(M, N); // ordinal rank matrix
    p.col(0) = Vec::LinSpaced(N, 0, N-1);
    r(0, p.col(0)) = Vec::LinSpaced(N, 0, N-1);

    for (auto i = 1; i < M; ++i) {
        p.col(i) = p.col(i - 1); // this is a critical part of the approach
        std::stable_sort(p.col(i).begin(), p.col(i).end(), [&](auto a, auto b) { return detail::less_than_f(points[a][i], points[b][i]); });
        r(i, p.col(i)) = Vec::LinSpaced(N, 0, N-1);
    }

    // 2) save min and max positions as well as the column index for the max position
    Vec maxc(N);
    Vec maxp(N);
    for (auto i = 0; i < N; ++i) {
        auto c = r.col(i);
        auto max = std::max_element(c.begin(), c.end());
        maxp(i) = *max;
        maxc(i) = std::distance(c.begin(), max);
    }
    // 3) compute ranks / fronts
    for (auto i : p(Eigen::seq(0, N - 2), 0)) {
        if (maxp(i) == N - 1) {
            continue;
        }

        for (auto j : p(Eigen::seq(maxp(i) + 1, N - 1), maxc(i))) {
            non_dom_rank[j] += static_cast<int64_t>(non_dom_rank[i] == non_dom_rank[j] && (r.col(i) < r.col(j)).all());
        }
    }
    auto n_fronts = *std::max_element(non_dom_rank.begin(), non_dom_rank.end()) + 1;
    non_dom_fronts.resize(n_fronts);
    for (auto i = 0; i < N; ++i) {
        non_dom_fronts[non_dom_rank[i]].push_back(i);
    }
    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}
} // namespace pagmo
