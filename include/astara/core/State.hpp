#ifndef ASTARA_CORE_STATE_HPP
#define ASTARA_CORE_STATE_HPP

/**
 * @file   State.hpp
 * @brief  Fixed- and dynamic-size state vectors with vector-space arithmetic.
 *
 * Most components in ASTARA have a small, compile-time-known number of state
 * variables (a reactor with N=15-20, a pressurizer with N=2, etc.).  We use
 * `std::array<double, N>` wrapped in a thin type that supplies the arithmetic
 * operators the integrators require:
 *
 *   - `operator+(State, State) -> State`
 *   - `operator-(State, State) -> State`
 *   - `operator*(State, double) -> State`
 *   - `operator+=`, `operator-=`, `operator*=`
 *
 * For dynamically-sized vectors (e.g. the integrator working over a flattened
 * concatenation of many components) the same operations are provided over
 * `std::vector<double>`.  Both share an interface so `Integrator<...>` can be
 * templated on the storage type.
 *
 * Why not Eigen?  Eigen would be the natural choice but pulling it as a hard
 * dependency for what amounts to "axpy on a 20-element vector" was judged
 * over-kill.  The operations below are inlinable and the compiler can
 * vectorise them as well as Eigen would.
 */

#include <array>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace astara::core {

// =============================================================================
// FixedState<N> -- compile-time-sized state vector
// =============================================================================

/**
 * @brief A compile-time-sized state vector with vector-space operators.
 * @tparam N  number of scalar state variables.
 *
 * `FixedState<N>` is a thin wrapper over `std::array<double, N>` that adds the
 * arithmetic the time integrators need.  It is trivially copyable, has no
 * heap allocation, and is `constexpr`-friendly where the standard library
 * permits.
 */
template <std::size_t N>
struct FixedState {
    std::array<double, N> data{};  ///< raw storage, zero-initialised

    /// Number of scalar state variables.
    static constexpr std::size_t size() noexcept { return N; }

    constexpr double  operator[](std::size_t i) const noexcept { return data[i]; }
    constexpr double& operator[](std::size_t i)       noexcept { return data[i]; }

    /// Iterators (so range-for and STL algorithms work).
    constexpr auto begin()       noexcept { return data.begin(); }
    constexpr auto end()         noexcept { return data.end();   }
    constexpr auto begin() const noexcept { return data.begin(); }
    constexpr auto end()   const noexcept { return data.end();   }

    /// Element-wise in-place addition.
    FixedState& operator+=(const FixedState& rhs) noexcept {
        for (std::size_t i = 0; i < N; ++i) data[i] += rhs.data[i];
        return *this;
    }
    /// Element-wise in-place subtraction.
    FixedState& operator-=(const FixedState& rhs) noexcept {
        for (std::size_t i = 0; i < N; ++i) data[i] -= rhs.data[i];
        return *this;
    }
    /// Scalar in-place multiplication.
    FixedState& operator*=(double s) noexcept {
        for (std::size_t i = 0; i < N; ++i) data[i] *= s;
        return *this;
    }

    /// Returns true iff every element is finite (not NaN, not inf).
    bool isFinite() const noexcept {
        for (std::size_t i = 0; i < N; ++i) {
            if (!std::isfinite(data[i])) return false;
        }
        return true;
    }
};

template <std::size_t N>
inline FixedState<N> operator+(FixedState<N> a, const FixedState<N>& b) noexcept { a += b; return a; }
template <std::size_t N>
inline FixedState<N> operator-(FixedState<N> a, const FixedState<N>& b) noexcept { a -= b; return a; }
template <std::size_t N>
inline FixedState<N> operator*(FixedState<N> a, double s) noexcept { a *= s; return a; }
template <std::size_t N>
inline FixedState<N> operator*(double s, FixedState<N> a) noexcept { a *= s; return a; }

// =============================================================================
// DynamicState -- dynamically-sized state vector
// =============================================================================

/**
 * @brief A run-time-sized state vector with vector-space operators.
 *
 * Used by the loop integrator when the number of state variables depends on
 * configuration (number of moderator nodes, number of SG nodes, etc.).
 * Otherwise prefer `FixedState<N>`.
 */
class DynamicState {
public:
    DynamicState() = default;
    explicit DynamicState(std::size_t n)            : data_(n, 0.0) {}
    DynamicState(std::size_t n, double v)           : data_(n, v)   {}
    explicit DynamicState(std::vector<double> data) : data_(std::move(data)) {}

    std::size_t size() const noexcept { return data_.size(); }
    void resize(std::size_t n)         { data_.resize(n, 0.0); }

    double  operator[](std::size_t i) const noexcept { assert(i < data_.size()); return data_[i]; }
    double& operator[](std::size_t i)       noexcept { assert(i < data_.size()); return data_[i]; }

    auto begin()       noexcept { return data_.begin(); }
    auto end()         noexcept { return data_.end();   }
    auto begin() const noexcept { return data_.begin(); }
    auto end()   const noexcept { return data_.end();   }

    const std::vector<double>& raw() const noexcept { return data_; }
    std::vector<double>&       raw()       noexcept { return data_; }

    DynamicState& operator+=(const DynamicState& rhs) {
        if (rhs.data_.size() != data_.size()) {
            throw std::invalid_argument("DynamicState size mismatch in +=");
        }
        for (std::size_t i = 0; i < data_.size(); ++i) data_[i] += rhs.data_[i];
        return *this;
    }
    DynamicState& operator-=(const DynamicState& rhs) {
        if (rhs.data_.size() != data_.size()) {
            throw std::invalid_argument("DynamicState size mismatch in -=");
        }
        for (std::size_t i = 0; i < data_.size(); ++i) data_[i] -= rhs.data_[i];
        return *this;
    }
    DynamicState& operator*=(double s) noexcept {
        for (auto& v : data_) v *= s;
        return *this;
    }

    bool isFinite() const noexcept {
        for (auto v : data_) if (!std::isfinite(v)) return false;
        return true;
    }

private:
    std::vector<double> data_;
};

inline DynamicState operator+(DynamicState a, const DynamicState& b) { a += b; return a; }
inline DynamicState operator-(DynamicState a, const DynamicState& b) { a -= b; return a; }
inline DynamicState operator*(DynamicState a, double s) noexcept   { a *= s; return a; }
inline DynamicState operator*(double s, DynamicState a) noexcept   { a *= s; return a; }

// =============================================================================
// Numerical helpers
// =============================================================================

/// Element-wise infinity-norm (max |a_i - b_i|), for tolerance comparisons.
template <std::size_t N>
inline double maxAbsDiff(const FixedState<N>& a, const FixedState<N>& b) noexcept {
    double m = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double d = std::abs(a[i] - b[i]);
        if (d > m) m = d;
    }
    return m;
}

inline double maxAbsDiff(const DynamicState& a, const DynamicState& b) {
    if (a.size() != b.size()) throw std::invalid_argument("size mismatch in maxAbsDiff");
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double d = std::abs(a[i] - b[i]);
        if (d > m) m = d;
    }
    return m;
}

/// A small absolute tolerance used in clamps and "is positive" checks.
inline constexpr double kTinyValue = 1e-12;

}  // namespace astara::core

#endif  // ASTARA_CORE_STATE_HPP
