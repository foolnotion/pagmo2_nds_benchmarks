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

#include <cpp-sort/sorters.h>

namespace pagmo
{

namespace detail {
    template<typename T>
    struct item {
        T value;
        size_t index;

        bool operator<(item other) const { return detail::less_than_f<T>(value, other.value); }
    };

    template<typename T /* block type */, size_t S = std::numeric_limits<T>::digits /* block size in bits */>
    class Bitset {
        std::vector<T> blocks_;
        size_t numBits_{};

        [[nodiscard]] inline auto BlockIndex(size_t i) const { return i / S; }
        [[nodiscard]] inline auto BitIndex(size_t i) const { return BlockIndex(i) + i % S; }

        public:
        static constexpr T ZeroBlock = T{0}; // a block with all the bits set to zero
        static constexpr T OneBlock = ~ZeroBlock; // a block with all the bits set to one
        static constexpr size_t BlockSize = S;
        using Block = T;

        Bitset() = default;

        explicit Bitset(size_t n, T blockInit)
        {
            Resize(n, blockInit);
        }

        inline void Fill(T value) { std::fill(blocks_.begin(), blocks_.end(), value); }

        inline void Resize(size_t n, T blockInit = T{0}) {
            numBits_ = n;
            size_t const nb = (n / S) + (n % S != 0); // NOLINT
            blocks_.resize(nb, blockInit);
            blocks_.back() >>= S * nb - n; // zero the bits in the last block that are over n
        }

        inline void Set(size_t i) {
            assert(i < numBits_);
            blocks_[i / S] |= (T{1} << (i % S));
        }

        inline void Reset(size_t i) {
            assert(i < numBits_);
            blocks_[i / S] &= ~(T{1} << (i % S));
        }

        [[nodiscard]] auto PopCount() const -> size_t {
            return std::transform_reduce(blocks_.begin(), blocks_.end(), size_t{0}, std::plus<>{}, [](auto b) { return __builtin_popcountl(b); });
        }

        auto Data() -> T* { return blocks_.data(); }
        [[nodiscard]] auto Data() const -> T const* { return blocks_.data(); }

        [[nodiscard]] inline auto NumBlocks() const -> size_t { return blocks_.size(); }

        [[nodiscard]] inline auto Empty() const -> bool { return blocks_.empty(); }

        inline auto operator[](size_t i) const -> bool {
            return static_cast<bool>(blocks_[i / S] & (T{1} << (i % S)));
        }

        friend auto operator&(Bitset const& lhs, Bitset const& rhs) -> Bitset {
            auto result = lhs;
            result &= rhs;
            return result;
        }

        friend auto operator&=(Bitset& lhs, Bitset const& rhs) -> Bitset& {
            assert(lhs.blocks_.size() == rhs.blocks_.size());
            assert(lhs.numBits_ == rhs.numBits_);
            for (size_t i = 0; i < lhs.blocks_.size(); ++i) {
                lhs[i] &= rhs[i];
            }
            return lhs;
        };

        friend auto operator<<(std::ostream& os, Bitset const& bs) -> std::ostream& {
            os << "{ ";
            for (auto i = 0UL; i < bs.blocks_.size(); ++i) {
                auto b = bs.blocks_[i];
                auto o = BlockSize * i;

                while (b) {
                    auto x = o + __builtin_ctzl(b);
                    b &= (b-1);
                    os << x;
                    if (b || (i < bs.blocks_.size() - 1 && bs.blocks_[i+1] != 0)) {
                        os << ", ";
                    }
                }
            }
            os << " }";
            return os;
        };

        [[nodiscard]] auto ToVec() const -> std::vector<size_t> {
            std::vector<size_t> result;
            result.reserve(PopCount());

            for (auto i = 0UL; i < blocks_.size(); ++i) {
                auto b = blocks_[i];
                auto o = BlockSize * i;

                while (b) {
                    auto x = o + __builtin_ctzl(b);
                    result.push_back(x);
                    b &= (b-1);
                }
            }

            return result;
        }

        template <typename U, std::enable_if_t<std::is_integral_v<U> && std::is_unsigned_v<U>, bool> = true>
        static inline auto CountTrailingZeros(U block) -> size_t // NOLINU
        {
            assert(block != U{0}); // output is undefined for 0!

            constexpr size_t u_bits_number = std::numeric_limits<unsigned>::digits;
            constexpr size_t ul_bits_number = std::numeric_limits<unsigned long>::digits;
            constexpr size_t ull_bits_number = std::numeric_limits<unsigned long long>::digits;
            const size_t bits_per_block = sizeof(U) * CHAR_BIT;

#if defined(__clang__) || defined(__GNUC__)
            if constexpr (bits_per_block <= u_bits_number) {
                return static_cast<size_t>(__builtin_ctz(static_cast<unsigned int>(block)));
            } else if constexpr (bits_per_block <= ul_bits_number) {
                return static_cast<size_t>(__builtin_ctzl(static_cast<unsigned long>(block)));
            } else if constexpr (bits_per_block <= ull_bits_number) {
                return static_cast<size_t>(__builtin_ctzll(static_cast<unsigned long long>(block)));
            }
#elif defined(_MSC_VER)
            constexpr size_t ui64_bits_number = std::numeric_limits<unsigned __int64>::digits;
            if constexpr (bits_per_block <= ul_bits_number) {
                unsigned long index = std::numeric_limits<unsigned long>::max();
                _BitScanForward(&index, static_cast<unsigned long>(block));
                return static_cast<size_t>(index);
            } else if constexpr (bits_per_block <= ui64_bits_number) {
#if (defined(_M_X64) || defined(_M_ARM64))
                unsigned long index = std::numeric_limits<unsigned long>::max();
                _BitScanForward64(&index, static_cast<unsigned __int64>(block));
                return static_cast<size_t>(index);
#else
                constexpr unsigned long max_ul = std::numeric_limits<unsigned long>::max();
                unsigned long low = block & max_ul;
                if (low != 0) {
                    unsigned long index = std::numeric_limits<unsigned long>::max();
                    _BitScanForward(&index, low);
                    return static_cast<size_t>(index);
                }
                unsigned long high = block >> ul_bits_number;
                unsigned long index = std::numeric_limits<unsigned long>::max();
                _BitScanForward(&index, high);
                return static_cast<size_t>(ul_bits_number + index);
#endif
            }
#endif
            U mask = U { 1 };
            for (size_t i = 0; i < bits_per_block; ++i) {
                if ((block & mask) != 0) {
                    return i;
                }
                mask <<= 1;
            }
        }
    };
}

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

    using Bitset = detail::Bitset<uint64_t>;

    Bitset b(N, Bitset::OneBlock);
    std::vector<Bitset> bs(N);                          // vector of bitsets (one for each individual)
    std::vector<std::pair<size_t, size_t>> br(N);       // vector of ranges keeping track of the first/last non-zero blocks
    std::vector<detail::item<vector_double::value_type>> items(N); // these items hold the values to be sorted along with their associated index

    auto const numBlocks = b.NumBlocks();
    auto const blockSize = Bitset::BlockSize;
    for (size_t i = 0; i < N; ++i) {
        items[i].index = i;
        b.Reset(i);
        bs[i] = b;
        br[i] = {0, numBlocks - 1};
    }
    std::vector<Bitset> rk;         // vector of sets keeping track of individuals whose rank was updated
    rk.emplace_back(N, Bitset::OneBlock);

    std::vector<size_t> rank(N, 0);

    cppsort::merge_sorter sorter;

    for (size_t k = 1; k < M; ++k) {
        for (size_t i = 0; i < N; ++i) {
            auto& item = items[i];
            item.value = points[item.index][k];
        }
        sorter(items);
        //std::stable_sort(items.begin(), items.end());
        b.Fill(Bitset::OneBlock);

        for (auto [_, i] : items) {
            b.Reset(i);
            auto [lo, hi] = br[i];
            if (lo > hi) { continue; }

            auto* p = bs[i].Data();
            auto const* q = b.Data();

            // tighten the interval around non-zero blocks
            while(lo <= hi && !(p[lo] & q[lo])) { ++lo; } // NOLINT
            while(lo <= hi && !(p[hi] & q[hi])) { --hi; } // NOLINT
            br[i] = {lo, hi};

            if (k < M-1) {
                // perform the set intersection
                for (size_t j = lo; j <= hi; ++j) {
                    p[j] &= q[j];
                }
            } else {
                auto rnk = non_dom_rank[i];
                if (rnk + 1UL == rk.size()) {
                    rk.emplace_back(N, Bitset::ZeroBlock);
                }
                auto* r = rk[rnk].Data();
                auto* s = rk[rnk + 1].Data();

                for (size_t j = lo; j <= hi; ++j) {
                    auto v = p[j] & q[j] & r[j]; // obtain the dominance set
                    r[j] &= ~v;                  // remove dominated individuals from current rank set
                    s[j] |= v;                   // add the individuals to the next rank set

                    auto o = blockSize * j;
                    while(v) {
                        auto x = o + Bitset::CountTrailingZeros(v);
                        v &= (v-1);
                        ++non_dom_rank[x];
                    }
                }
            }
        }
    }

    auto n_fronts = *std::max_element(non_dom_rank.begin(), non_dom_rank.end()) + 1;
    non_dom_fronts.resize(n_fronts);
    for (size_t i = 0UL; i < N; ++i) {
        non_dom_fronts[non_dom_rank[i]].push_back(i);
    }
    return std::make_tuple(std::move(non_dom_fronts), std::move(dom_list), std::move(dom_count),
                           std::move(non_dom_rank));
}
} // namespace pagmo
