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

#ifndef PAGMO_CEC2013_DATA_HPP
#define PAGMO_CEC2013_DATA_HPP

#include <unordered_map>
#include <vector>

namespace pagmo
{

namespace detail
{

namespace cec2013_data
{

extern const std::unordered_map<unsigned, std::vector<double>> MD;
extern const std::vector<double> shift_data;

} // namespace cec2013_data

} // namespace detail

} // namespace pagmo

#endif
