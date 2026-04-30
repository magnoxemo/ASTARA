/**
 * @file   test_state.cpp
 * @brief  Unit tests for FixedState<N> and DynamicState arithmetic.
 */

#include "astara/core/State.hpp"
#include <gtest/gtest.h>

using astara::core::FixedState;
using astara::core::DynamicState;

TEST(FixedState, DefaultsToZero) {
    FixedState<4> s;
    for (std::size_t i = 0; i < s.size(); ++i) EXPECT_EQ(s[i], 0.0);
    EXPECT_TRUE(s.isFinite());
}

TEST(FixedState, AdditionAndSubtraction) {
    FixedState<3> a{{1.0, 2.0, 3.0}};
    FixedState<3> b{{4.0, 5.0, 6.0}};
    auto c = a + b;
    EXPECT_DOUBLE_EQ(c[0], 5.0);
    EXPECT_DOUBLE_EQ(c[1], 7.0);
    EXPECT_DOUBLE_EQ(c[2], 9.0);
    auto d = c - b;
    EXPECT_DOUBLE_EQ(d[0], 1.0);
    EXPECT_DOUBLE_EQ(d[1], 2.0);
    EXPECT_DOUBLE_EQ(d[2], 3.0);
}

TEST(FixedState, ScalarMultiplicationCommutes) {
    FixedState<2> a{{2.0, -3.0}};
    auto b = a * 4.0;
    auto c = 4.0 * a;
    EXPECT_DOUBLE_EQ(b[0], c[0]);
    EXPECT_DOUBLE_EQ(b[1], c[1]);
    EXPECT_DOUBLE_EQ(b[0], 8.0);
    EXPECT_DOUBLE_EQ(b[1], -12.0);
}

TEST(FixedState, IsFiniteDetectsNaN) {
    FixedState<3> s{{1.0, std::nan(""), 3.0}};
    EXPECT_FALSE(s.isFinite());
}

TEST(DynamicState, ResizeAndArithmetic) {
    DynamicState a(4, 1.5);
    DynamicState b(4, 0.5);
    auto c = a + b;
    EXPECT_EQ(c.size(), 4u);
    for (std::size_t i = 0; i < 4; ++i) EXPECT_DOUBLE_EQ(c[i], 2.0);

    a *= 2.0;
    for (std::size_t i = 0; i < 4; ++i) EXPECT_DOUBLE_EQ(a[i], 3.0);
}

TEST(DynamicState, SizeMismatchThrows) {
    DynamicState a(3, 0.0);
    DynamicState b(4, 0.0);
    EXPECT_THROW({ a += b; }, std::invalid_argument);
}
