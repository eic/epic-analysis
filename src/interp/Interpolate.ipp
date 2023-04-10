// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#ifndef Interpolate_IPP_
#define Interpolate_IPP_

#include "Interpolate.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

// Version of `modf` with correct sign behaviour.
template<typename T>
T modf(T x, T* iptr) {
	if (!(x < 0.)) {
		return std::modf(x, iptr);
	} else {
		T frac = std::modf(x, iptr);
		frac += 1;
		*iptr -= 1;
		return frac;
	}
}

template<typename T, std::size_t N>
GridView<T, N>::GridView(T const* data, Coords<T, N> coords) :
		_data(data),
		_coords(coords) {
	_count_total = 1;
	for (std::size_t idx = 0; idx < N; ++idx) {
		_count_total *= _coords[idx].size();
		if (_coords[idx].size() <= 1) {
			throw SingularDimensionError(idx);
		}
	}
}

template<typename T, std::size_t N>
GridView<T, N - 1> GridView<T, N>::operator[](std::size_t idx) const {
	std::size_t width = _count_total / count(0);
	std::array<std::vector<T>, N - 1> coords_red;
	std::copy(_coords.begin() + 1, _coords.end(), coords_red.begin());
	return GridView<T, N - 1>(_data + idx * width, coords_red);
}

template<typename T, std::size_t N>
T const& GridView<T, N>::operator[](CellIndex<N> idx) const {
	CellIndex<N - 1> idx_red;
	std::copy(idx.begin() + 1, idx.end(), idx_red.begin());
	return (*this)[idx[0]][idx_red];
}

template<typename T, std::size_t N>
Grid<T, N>::Grid(T const* data, Coords<T, N> coords) :
		_data(),
		_coords(coords) {
	_count_total = 1;
	for (std::size_t idx = 0; idx < N; ++idx) {
		_count_total *= _coords[idx].size();
		if (_coords[idx].size() <= 1) {
			throw SingularDimensionError(idx);
		}
	}
	_data = std::vector<T>(data, data + _count_total);
}

template<typename T, std::size_t N>
T GridView<T, N>::space_from_grid(std::size_t dim, std::size_t idx) const {
	return _coords[dim][idx];
}
template<typename T, std::size_t N>
std::size_t GridView<T, N>::grid_from_space(std::size_t dim, T x) const {
	T frac;
	return grid_from_space(dim, x, &frac);
}
template<typename T, std::size_t N>
std::size_t GridView<T, N>::grid_from_space(std::size_t dim, T x, T* frac) const {
	std::size_t lower = 0;
	std::size_t upper = _coords[dim].size() - 1;
	T lower_x = _coords[dim][lower];
	T upper_x = _coords[dim][upper];
	while (lower + 1 < upper) {
		T diff = (x - lower_x) / (upper_x - lower_x);
		if (!(diff >= 0.) || !(diff <= 1.)) {
			return std::numeric_limits<std::size_t>::max();
		}
		T idx_float;
		modf(diff * (upper - lower - 1), &idx_float);
		std::size_t idx = static_cast<std::size_t>(idx_float) + lower;
		T lower_x_new = _coords[dim][idx];
		T upper_x_new = _coords[dim][idx + 1];
        if (x >= lower_x_new && x < upper_x_new) {
            lower = idx;
            upper = idx + 1;
            lower_x = lower_x_new;
            upper_x = upper_x_new;
        } else if (x >= upper_x_new) {
            lower = idx + 1;
            lower_x = upper_x_new;
        } else if (x < lower_x_new) {
            upper = idx;
            upper_x = lower_x_new;
        }
	}
	*frac = x - lower_x;
	return lower;
}

template<typename T, std::size_t N>
Point<T, N> GridView<T, N>::space_from_grid(CellIndex<N> idx) const {
	Point<T, N> result;
	for (std::size_t dim = 0; dim < N; ++dim) {
		result[dim] = space_from_grid(dim, idx[dim]);
	}
	return result;
}
template<typename T, std::size_t N>
CellIndex<N> GridView<T, N>::grid_from_space(Point<T, N> x) const {
	CellIndex<N> result;
	for (std::size_t dim = 0; dim < N; ++dim) {
		result[dim] = grid_from_space(dim, x[dim]);
	}
	return result;
}
template<typename T, std::size_t N>
CellIndex<N> GridView<T, N>::grid_from_space(Point<T, N> x, Point<T, N>* frac) const {
	CellIndex<N> result;
	for (std::size_t dim = 0; dim < N; ++dim) {
		result[dim] = grid_from_space(dim, x[dim], &frac[dim]);
	}
	return result;
}

template<typename T, std::size_t N>
T LinearView<T, N>::operator()(Point<T, N> x) const {
	T frac;
	std::size_t idx = _grid.grid_from_space(0, x[0], &frac);
	Point<T, N - 1> x_next;
	std::copy(x.begin() + 1, x.end(), x_next.begin());
	if (idx >= 0 && idx < _grid.count(0) - 1) {
		T linear_1 = LinearView<T, N - 1>(_grid.reduce(idx))(x_next);
		T linear_2 = LinearView<T, N - 1>(_grid.reduce(idx + 1))(x_next);
		T delta = _grid.space_from_grid(0, idx + 1) - _grid.space_from_grid(0, idx);
		return linear(linear_1, linear_2, delta, frac);
	} else {
		return std::numeric_limits<T>::quiet_NaN();
	}
}

template<typename T, std::size_t N>
T CubicView<T, N>::operator()(Point<T, N> x) const {
	T frac;
	std::size_t idx = _grid.grid_from_space(0, x[0], &frac);
	Point<T, N - 1> x_next;
	std::copy(x.begin() + 1, x.end(), x_next.begin());
	if (idx >= 0 && idx < _grid.count(0) - 1) {
		T cubic_0 = idx == 0 ? 0. :
			CubicView<T, N - 1>(_grid.reduce(idx - 1))(x_next);
		T cubic_1 = CubicView<T, N - 1>(_grid.reduce(idx))(x_next);
		T cubic_2 = CubicView<T, N - 1>(_grid.reduce(idx + 1))(x_next);
		T cubic_3 = idx >= _grid.count(0) - 2 ? 0. :
			CubicView<T, N - 1>(_grid.reduce(idx + 2))(x_next);
		T delta_0 = idx == 0 ? std::numeric_limits<T>::infinity() :
			_grid.space_from_grid(0, idx) - _grid.space_from_grid(0, idx - 1);
		T delta_1 = _grid.space_from_grid(0, idx + 1) - _grid.space_from_grid(0, idx);
		T delta_2 = idx >= _grid.count(0) - 2 ? std::numeric_limits<T>::infinity() :
			_grid.space_from_grid(0, idx + 2) - _grid.space_from_grid(0, idx + 1);
		return cubic(cubic_0, cubic_1, cubic_2, cubic_3, delta_0, delta_1, delta_2, frac);
	} else {
		return std::numeric_limits<T>::quiet_NaN();
	}
}

template<typename T, std::size_t N, std::size_t K>
inline std::array<Grid<T, N>, K> read_grids(
		std::vector<std::array<T, N + K> > const& raw_data) {
	std::vector<Point<T, N> > grid_points(raw_data.size());
	std::vector<std::array<T, K> > data(raw_data.size());
	for (std::size_t row_idx = 0; row_idx < raw_data.size(); ++row_idx) {
		for (std::size_t dim_idx = 0; dim_idx < N; ++dim_idx) {
			grid_points[row_idx][dim_idx] = raw_data[row_idx][dim_idx];
		}
		for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
			data[row_idx][col_idx] = raw_data[row_idx][N + col_idx];
		}
	}

	// If there are not enough grid points to form one unit cell, then give up.
	if (grid_points.size() < (1 << N)) {
		throw NotEnoughPointsError(grid_points.size(), 1 << N);
	}

	std::size_t count_total = 1;
	std::array<bool, N> filled_cols = { false };
	std::array<std::size_t, N> dim_permute_map;
	CellIndex<N> count;
	std::array<std::vector<T>, N> coords;
	for (std::size_t dim_idx = 0; dim_idx < N; ++dim_idx) {
		std::size_t sub_count = grid_points.size() / count_total;
		if (sub_count * count_total != grid_points.size()) {
			throw NotEnoughPointsError(
				grid_points.size(),
				(sub_count + 1) * count_total);
		}
		// Find the column that steps by `count_total`.
		std::size_t step_dim_idx = 0;
		while (true) {
			if (step_dim_idx == N) {
				throw NotEnoughPointsError(
					grid_points.size(),
					2 * grid_points.size());
			}
			if (filled_cols[step_dim_idx]) {
				step_dim_idx += 1;
				continue;
			}
			T initial_point = grid_points[0][step_dim_idx];
			T final_point = grid_points[count_total][step_dim_idx];
			if (initial_point == final_point) {
				step_dim_idx += 1;
				continue;
			}
			bool same = true;
			for (std::size_t row_idx = 0; row_idx < count_total; ++row_idx) {
				if (grid_points[row_idx][step_dim_idx] != initial_point) {
					same = false;
					break;
				}
			}
			if (!same) {
				step_dim_idx += 1;
				continue;
			}
			filled_cols[step_dim_idx] = true;
			dim_permute_map[dim_idx] = step_dim_idx;
			break;
		}

		// Use the column to fill a vector of the grid spacings.
		std::vector<T> next_coords;
		for (std::size_t group_idx = 0; group_idx < sub_count; ++group_idx) {
			T next = grid_points[group_idx * count_total][step_dim_idx];
			if (next_coords.empty() || next != next_coords[0]) {
				next_coords.push_back(next);
			} else {
				break;
			}
		}
		std::size_t next_count = next_coords.size();
		count[step_dim_idx] = next_count;
		if (next_count <= 1) {
			throw SingularDimensionError(step_dim_idx);
		}

		// Check that the data in the column is valid.
		for (std::size_t row_idx = 0; row_idx < grid_points.size(); ++row_idx) {
			T next = grid_points[row_idx][step_dim_idx];
			std::size_t group_idx = row_idx / count_total;
			if (next != next_coords[group_idx % next_count]) {
				throw UnexpectedGridPointError(row_idx);
			}
		}

		// Check that the grid is increasing.
		for (std::size_t idx = 1; idx < next_coords.size(); ++idx) {
			if (!(next_coords[idx] > next_coords[idx - 1])) {
				throw InvalidGridPlanesError();
			}
		}
		coords[step_dim_idx] = next_coords;

		// Increase count for next round.
		count_total *= next_count;
	}

	// Use the grid information to reorder the data.
	std::array<std::vector<T>, K> data_transposed;
	std::array<std::size_t, N> count_totals = { count_total / count[0] };
	for (std::size_t dim_idx = 1; dim_idx < N; ++dim_idx) {
		count_totals[dim_idx] = count_totals[dim_idx - 1] / count[dim_idx];
	}
	for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
		data_transposed[col_idx].resize(data.size());
	}
	for (std::size_t row_idx = 0; row_idx < data.size(); ++row_idx) {
		std::size_t new_row_idx = 0;
		std::size_t new_count = 1;
		for (std::size_t dim_idx = 0; dim_idx < N; ++dim_idx) {
			std::size_t new_dim_idx = dim_permute_map[dim_idx];
			std::size_t rel_idx = (row_idx / new_count) % count[new_dim_idx];
			std::size_t old_count = count_totals[new_dim_idx];
			new_row_idx += rel_idx * old_count;
			new_count *= count[new_dim_idx];
		}
		for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
			data_transposed[col_idx][new_row_idx] = data[row_idx][col_idx];
		}
	}

	std::array<Grid<T, N>, K> grids;
	for (std::size_t col_idx = 0; col_idx < K; ++col_idx) {
		grids[col_idx] = Grid<T, N>(data_transposed[col_idx].data(), coords);
	}
	return grids;
}

inline NotEnoughPointsError::NotEnoughPointsError(
	std::size_t points,
	std::size_t expected_points) :
	std::runtime_error(
		"Not enough points to construct the grid (found "
		+ std::to_string(points) + ", expected "
		+ std::to_string(expected_points) + ")"),
	points(points),
	expected_points(expected_points) { }

inline SingularDimensionError::SingularDimensionError(std::size_t dim) :
	std::runtime_error("Grid is singular in dimension " + std::to_string(dim)),
	dim(dim) { }

inline InvalidGridPlanesError::InvalidGridPlanesError() :
	std::runtime_error("Lower grid bound must be smaller than upper bound") { }

inline InvalidSpacingError::InvalidSpacingError() :
	std::runtime_error("Grid must be spaced uniformly") { }

inline UnexpectedGridPointError::UnexpectedGridPointError(
	std::size_t line_number) :
	std::runtime_error(
		"Unexpected grid point at line " + std::to_string(line_number)),
	line_number(line_number) { }

#endif

