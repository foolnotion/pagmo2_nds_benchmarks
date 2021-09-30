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

#ifndef PAGMO_NONDOMINATED_SORTING_HPP
#define PAGMO_NONDOMINATED_SORTING_HPP

#include <cstddef>
#include <climits>
#include <limits>

#include <pagmo/detail/visibility.hpp>
#include <pagmo/exceptions.hpp>
#include <pagmo/types.hpp>

namespace pagmo
{
    enum class non_dominated_sorting_algorithm {
        DS,     // deductive sort
        ENS_BS, // binary efficient non-dominated sort
        ENS_SS, // sequential efficient non-dominated sort
        HS,     // hierarchical sort
        MNDS,   // merge non-dominated sort
        NDS,    // fast non-dominated sort
        RS      // rank sort
    };

namespace detail
{
    template <typename T, std::enable_if_t<std::is_integral_v<T> && std::is_unsigned_v<T>, bool> = true>
    static inline size_t count_trailing_zeros(T block)
    {
        assert(block != T{0}); // output is undefined for 0!

        constexpr size_t u_bits_number = std::numeric_limits<unsigned>::digits;
        constexpr size_t ul_bits_number = std::numeric_limits<unsigned long>::digits;
        constexpr size_t ull_bits_number = std::numeric_limits<unsigned long long>::digits;
        const size_t bits_per_block = sizeof(T) * CHAR_BIT;

#if defined(__clang__) || defined(__GNUC__)
        if constexpr (bits_per_block <= u_bits_number) {
            return static_cast<size_t>(__builtin_ctz(static_cast<unsigned int>(block)));
        } else if constexpr (bits_per_block <= ul_bits_number) {
            return static_cast<size_t>(__builtin_ctzl(static_cast<unsigned long>(block)));
        } else if constexpr (bits_per_block <= ull_bits_number) {
            return static_cast<size_t>(__builtin_ctzll(static_cast<unsigned long long>(block)));
        }
#elif defined(_MSC_VER)
#if defined(_M_IX86) || defined(_M_ARM) || defined(_M_X64) || defined(_M_ARM64)
#include <intrin.h>
#pragma intrinsic(_BitScanForward)
        constexpr size_t ul_bits_number = std::numeric_limits<unsigned long>::digits;
        constexpr size_t ui64_bits_number = std::numeric_limits<unsigned __int64>::digits;
        if constexpr (bits_per_block <= ul_bits_number) {
            unsigned long index = std::numeric_limits<unsigned long>::max();
            _BitScanForward(&index, static_cast<unsigned long>(block));
            return static_cast<size_type>(index);
        } else if constexpr (bits_per_block <= ui64_bits_number) {
#if (defined(_M_X64) || defined(_M_ARM64))
#pragma intrinsic(_BitScanForward64)
            unsigned long index = std::numeric_limits<unsigned long>::max();
            _BitScanForward64(&index, static_cast<unsigned __int64>(block));
            return static_cast<size_type>(index);
#else
            constexpr unsigned long max_ul = std::numeric_limits<unsigned long>::max();
            unsigned long low = block & max_ul;
            if (low != 0) {
                unsigned long index = std::numeric_limits<unsigned long>::max();
                _BitScanForward(&index, low);
                return static_cast<size_type>(index);
            }
            unsigned long high = block >> ul_bits_number;
            unsigned long index = std::numeric_limits<unsigned long>::max();
            _BitScanForward(&index, high);
            return static_cast<size_type>(ul_bits_number + index);
#endif
        }
#endif
#endif
        T mask = T { 1 };
        for (size_t i = 0; i < bits_per_block; ++i) {
            if ((block & mask) != 0) {
                return i;
            }
            mask <<= 1;
        }
    }
} // namespace detail

/// Return type for the fast_non_dominated_sorting algorithm
using fnds_return_type = std::tuple<std::vector<std::vector<pop_size_t>>, std::vector<std::vector<pop_size_t>>,
      std::vector<pop_size_t>, std::vector<pop_size_t>>;

// Wrapper method which is parameterized with the concrete sorting algorithm to use and a duplicates-handling strategy
PAGMO_DLL_PUBLIC fnds_return_type non_dominated_sorting(const std::vector<vector_double> &, non_dominated_sorting_algorithm, bool = false /* dominate on equal */);

// Pareto-dominance
PAGMO_DLL_PUBLIC bool pareto_dominance(const vector_double &, const vector_double &);

// Non dominated front 2D (Kung's algorithm)
PAGMO_DLL_PUBLIC std::vector<pop_size_t> non_dominated_front_2d(const std::vector<vector_double> &);

// Fast non dominated sorting
PAGMO_DLL_PUBLIC fnds_return_type fast_non_dominated_sorting(const std::vector<vector_double> &);

// Deductive sorting
PAGMO_DLL_PUBLIC fnds_return_type deductive_sorting(const std::vector<vector_double> &);

// Hierarchical sorting
PAGMO_DLL_PUBLIC fnds_return_type hierarchical_sorting(const std::vector<vector_double> &);

// Efficient sorting (binary) 
PAGMO_DLL_PUBLIC fnds_return_type efficient_sorting_binary(const std::vector<vector_double> &);

// Efficient sorting (sequential) 
PAGMO_DLL_PUBLIC fnds_return_type efficient_sorting_sequential(const std::vector<vector_double> &);

// Merge non-dominated sorting
PAGMO_DLL_PUBLIC fnds_return_type merge_non_dominated_sorting(const std::vector<vector_double> &);

// Rank non-dominated sorting
PAGMO_DLL_PUBLIC fnds_return_type rank_non_dominated_sorting(const std::vector<vector_double> &);

} // namespace pagmo
#endif
