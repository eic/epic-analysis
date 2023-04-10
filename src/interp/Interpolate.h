// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#ifndef Interpolate_
#define Interpolate_

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

/// Linearly interpolate between two endpoints, with \p x between 0 and 1.
template<typename T>
inline T linear(T f1, T f2, T x) {
	return (1. - x) * f1 + x * f2;
}

/// Linear interpolation for non-uniformly spaced points.
template<typename T>
inline T linear(T f1, T f2, T h, T x) {
	return ((h - x) / h) * f1 + (x / h) * f2;
}

/// Cubic interpolation using \p f2 and \p f3 as the endpoints, and \p f1 and
/// \p f4 as the anchor points beyond the endpoints. \p x between 0 and 1.
template<typename T>
inline T cubic(T f0, T f1, T f2, T f3, T x) {
	return
		2. * (1. - x) * (1. - x) * (0.5 + x) * f1
		+ 0.5 * (1. - x) * (1. - x) * x * (f2 - f0)
		+ 2. * x * x * (1.5 - x) * f2
		- 0.5 * x * x * (1. - x) * (f3 - f1);
}

/// Cubic interpolation for non-uniformly spaced points.
template<typename T>
inline T cubic(T f0, T f1, T f2, T f3, T h0, T h1, T h2, T x) {
	T h01 = h0 + h1;
	T h12 = h1 + h2;
	T h = h0 + h1 + h2;
	T c0, c1, c2, c3;
	if (!std::isfinite(h0) && !std::isfinite(h2)) {
		return linear(f1, f2, h1, x);
	} else if (!std::isfinite(h0)) {
		c0 = 0.;
		c1 = (h1 - x) * (h1 * h12 - x * x) / (h1 * h1 * h12);
		c2 = x * (h1 * h2 + h1 * x - x * x) / (h2 * h1 * h1);
		c3 = -x * x * (h1 - x) / (h1 * h2 * h12);
	} else if (!std::isfinite(h2)) {
		c0 = -(h1 - x) * (h1 - x) * x / (h0 * h1 * h01);
		c1 = (h1 - x) * (h0 * h1 + h1 * x - x * x) / (h0 * h1 * h1);
		c2 = x * (h0 * h1 + 2. * h1 * x - x * x) / (h1 * h1 * h01);
		c3 = 0.;
	} else {
		c0 = -(h1 - x) * (h1 - x) * x / (h0 * h1 * h01);
		c1 = (h1 - x) * (h0 * h1 * h12 + h1 * h12 * x - h * x * x) / (h0 * h1 * h1 * h12);
		c2 = x * (h0 * h1 * h2 + h1 * (h + h2) * x - h * x * x) / (h2 * h1 * h1 * h01);
		c3 = -x * x * (h1 - x) / (h1 * h2 * h12);
	}
	return c0 * f0 + c1 * f1 + c2 * f2 + c3 * f3;
}

template<std::size_t N>
using CellIndex = std::array<std::size_t, N>;
template<typename T, std::size_t N>
using Point = std::array<T, N>;
template<typename T, std::size_t N>
using Coords = std::array<std::vector<T>, N>;

template<typename T, std::size_t N>
class Grid;

enum class Scale {
	LINEAR = 0,
	LOG    = 1,
};

/**
 * A view into an Grid (or another dataset). The difference between a GridView
 * and a Grid is that a Grid owns its data, while a GridView provides a Grid-
 * like interface into an already existing set of data.
 */
template<typename T, std::size_t N>
class GridView final {
	T const* _data;
	Coords<T, N> _coords;
	std::size_t _count_total;

public:
	/// Construct an GridView out of \p data. The data points are located at the
	/// positions given in the \p grid array.
	GridView(T const* data, Coords<T, N> coords);
	//GridView(T const* data, std::size_t count[N], T lower[N], T upper[N]);
	//GridView(T const* data, std::size_t count[N], T lower[N], T upper[N], Scale scale[N]);

	/// Number of data points in dimension \p dim.
	std::size_t count(std::size_t dim) const {
		return _coords[dim].size();
	}
	/// Number of data points in each dimension.
	CellIndex<N> count() const {
		CellIndex<N> result;
		for (std::size_t dim = 0; dim < N; ++dim) {
			result[dim] = _coords[dim].size();
		}
		return result;
	}
	/// Total number of data points.
	std::size_t count_total() const {
		return _count_total;
	}
	/// Access to the contigious underlying data.
	T const* data() const {
		return _data;
	}

	/// Access a slice of the data.
	GridView<T, N - 1> operator[](std::size_t idx) const;
	/// Access an element from the grid.
	T const& operator[](CellIndex<N> idx) const;
	/// Reduce the dimension of the GridView by accessing a slice.
	GridView<T, N - 1> reduce(std::size_t idx) const {
		return operator[](idx);
	}

	/// Get the grid index associated with a particular point in space.
	std::size_t grid_from_space(std::size_t dim, T x) const;
	std::size_t grid_from_space(std::size_t dim, T x, T* frac) const;
	CellIndex<N> grid_from_space(Point<T, N> x) const;
	CellIndex<N> grid_from_space(Point<T, N> x, Point<T, N>* frac) const;

	/// Gets the point in space associated with a grid index.
	T space_from_grid(std::size_t dim, std::size_t idx) const;
	Point<T, N> space_from_grid(CellIndex<N> idx) const;
};

/// Base specialization of GridView.
template<typename T>
class GridView<T, 0> {
	T const* _data;

public:
	explicit GridView(T const* data) : _data(data) { }
	GridView(T const* data, Coords<T, 0> coords) :
			_data(data) {
		static_cast<void>(coords);
	}

	CellIndex<0> count() const {
		return CellIndex<0>();
	}
	std::size_t count_total() const {
		return 1;
	}
	T const* data() const {
		return _data;
	}

	T const& operator[](CellIndex<0> idx) {
		return *_data;
	}
	operator T const&() const {
		return *_data;
	}
};

/**
 * A set of N-d data. A Grid is an owned version of a GridView.
 *
 * \sa GridView
 */
template<typename T, std::size_t N>
class Grid final {
	std::vector<T> _data;
	Coords<T, N> _coords;
	std::size_t _count_total;

public:
	Grid() : _data(), _coords(), _count_total(0) { }
	Grid(T const* data, Coords<T, N> coords);
	operator GridView<T, N>() const {
		return GridView<T, N>(_data.data(), _coords);
	}

	/// \copydoc GridView::count()
	std::size_t count(std::size_t dim) const {
		return static_cast<GridView<T, N> >(*this).count(dim);
	}
	/// \copydoc GridView::count()
	CellIndex<N> count() const {
		return static_cast<GridView<T, N> >(*this).count();
	}
	/// \copydoc GridView::count_total()
	std::size_t count_total() const {
		return _count_total;
	}
	/// \copydoc GridView::data()
	T const* data() const {
		return _data.data();
	}

	/// \copydoc GridView::operator[]()
	GridView<T, N - 1> operator[](std::size_t idx) const {
		return static_cast<GridView<T, N> >(*this)[idx];
	}
	/// \copydoc GridView::operator[]()
	T const& operator[](CellIndex<N> idx) const {
		return static_cast<GridView<T, N> >(*this)[idx];
	}
	/// \copydoc GridView::reduce()
	GridView<T, N - 1> reduce(std::size_t idx) const {
		return operator[](idx);
	}

	/// \copydoc GridView::grid_from_space()
	std::size_t grid_from_space(std::size_t dim, T x) const {
		return static_cast<GridView<T, N> >(*this).grid_from_space(dim, x);
	}
	std::size_t grid_from_space(std::size_t dim, T x, T* frac) const {
		return static_cast<GridView<T, N> >(*this).grid_from_space(dim, x, frac);
	}
	CellIndex<N> grid_from_space(Point<T, N> x) const {
		return static_cast<GridView<T, N> >(*this).grid_from_space(x);
	}
	CellIndex<N> grid_from_space(Point<T, N> x, Point<T, N>* frac) const {
		return static_cast<GridView<T, N> >(*this).grid_from_space(x, frac);
	}

	/// \copydoc GridView::space_from_grid()
	T space_from_grid(std::size_t dim, std::size_t idx) const {
		return static_cast<GridView<T, N> >(*this).space_from_grid(dim, idx);
	}
	Point<T, N> space_from_grid(CellIndex<N> idx) const {
		return static_cast<GridView<T, N> >(*this).space_from_grid(idx);
	}
};

/// Base specialization of Grid.
template<typename T>
class Grid<T, 0> final {
	T _data;

public:
	Grid() : _data() { }
	Grid(T const* data) : _data(*data) { };
	Grid(T const* data, Coords<T, 0> coords) : _data(*data) {
		static_cast<void>(coords);
	};
	operator GridView<T, 0>() const {
		return GridView<T, 0>(&_data);
	}

	/// \copydoc GridView::count()
	CellIndex<0> count() const {
		return CellIndex<0>();
	}
	/// \copydoc GridView::count_total()
	std::size_t count_total() const {
		return 1;
	}
	/// \copydoc GridView::data()
	T const* data() const {
		return &_data;
	}

	/// \copydoc GridView::operator[]()
	T const& operator[](CellIndex<0> idx) const {
		return static_cast<GridView<T, 0> >(*this)[idx];
	}
};

/**
 * A function that uses linear interpolation for a regularly spaced rectangular
 * grid of N-d data.
 */
template<typename T, std::size_t N>
class LinearView final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a LinearView of a non-floating-point type.");

	GridView<T, N> _grid;

public:
	LinearView(GridView<T, N> grid) : _grid(grid) { }
	/// Interpolate from the underlying GridView.
	T operator()(Point<T, N> x) const;
};

/// Base specialization of LinearView.
template<typename T>
class LinearView<T, 0> final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a LinearView of a non-floating-point type.");

	GridView<T, 0> _grid;

public:
	explicit LinearView(GridView<T, 0> grid) : _grid(grid) { }
	T operator()(Point<T, 0> x) const {
		static_cast<void>(x);
		return static_cast<T>(_grid);
	}
};

/**
 * A function that uses cubic interpolation for a regularly spaced rectangular
 * grid of N-d data.
 */
template<typename T, std::size_t N>
class CubicView final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a CubicView of a non-floating-point type.");

	GridView<T, N> _grid;

public:
	explicit CubicView(GridView<T, N> grid) : _grid(grid) { }
	/// Interpolate from the underlying GridView.
	T operator()(Point<T, N> x) const;
};

/// Base specialization of CubicView.
template<typename T>
class CubicView<T, 0> final {
	static_assert(
		std::is_floating_point<T>::value,
		"Cannot have a CubicView of a non-floating-point type.");

	GridView<T, 0> _grid;

public:
	explicit CubicView(GridView<T, 0> grid) : _grid(grid) { }
	T operator()(Point<T, 0> x) const {
		static_cast<void>(x);
		return static_cast<T>(_grid);
	}
};

/// Loads grids from an array of tuples. The provided data points must be
/// provided either in row-major order or column-major order. For each data
/// point, the first \p N numbers give the coordinates, and the next \p K
/// numbers give the Grid values at those coordinates. \p K different Grid%s
/// are returned.
template<typename T, std::size_t N, std::size_t K = 1>
std::array<Grid<T, N>, K> read_grids(
	std::vector<std::array<T, N + K> > const& raw_data);

/// Not enough points were provided to form a Grid without ragged edges.
struct NotEnoughPointsError : public std::runtime_error {
	std::size_t points;
	std::size_t expected_points;
	NotEnoughPointsError(std::size_t points, std::size_t expected_points);
};

/// One of the dimensions of a Grid is only one data point thick.
struct SingularDimensionError : public std::runtime_error {
	std::size_t dim;
	SingularDimensionError(std::size_t dim);
};

/// The hyper-cube bounds on a Grid invalid, most likely because the hyper-cube
/// has negative volume.
struct InvalidGridPlanesError : public std::runtime_error {
	InvalidGridPlanesError();
};

/// Data points used to construct a Grid are not spaced evenly.
struct InvalidSpacingError : public std::runtime_error {
	InvalidSpacingError();
};

/// Data points used to construct a Grid are provided in an unexpected order.
struct UnexpectedGridPointError : public std::runtime_error {
	std::size_t line_number;
	UnexpectedGridPointError(std::size_t line_number);
};

#include "Interpolate.ipp"

#endif

