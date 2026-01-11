#include "cmath"
#include <gtest/gtest.h>

#include "Parser.h"

TEST(AstaraParser, ConstantExpression) {
  astara::Function f("2 + 3 * 4", {});
  EXPECT_DOUBLE_EQ(f(), 14.0);
}

TEST(AstaraParser, Parentheses) {
  astara::Function f("(2 + 3) * 4", {});
  EXPECT_DOUBLE_EQ(f(), 20.0);
}

TEST(AstaraParser, DivisionAndSubtraction) {
  astara::Function f("10 - 6 / 2", {});
  EXPECT_DOUBLE_EQ(f(), 7.0);
}

TEST(AstaraParser, UnaryMinus) {
  astara::Function f("-5 + 2", {});
  EXPECT_DOUBLE_EQ(f(), -3.0);
}

TEST(AstaraParser, DoubleUnary) {
  astara::Function f("--5", {});
  EXPECT_DOUBLE_EQ(f(), 5.0);
}

TEST(AstaraParser, PowerOperator) {
  astara::Function f("2 ^ 3", {});
  EXPECT_DOUBLE_EQ(f(), 8.0);
}

TEST(AstaraParser, PowerRightAssociative) {
  astara::Function f("2 ^ 3 ^ 2", {});
  EXPECT_DOUBLE_EQ(f(), 512.0); // 2^(3^2)
}

TEST(AstaraParser, SingleVariable) {
  astara::Function f("x * 2", {"x"});
  EXPECT_DOUBLE_EQ(f(3.0), 6.0);
}

TEST(AstaraParser, MultipleVariables) {
  astara::Function f("x + y * 2", {"x", "y"});
  EXPECT_DOUBLE_EQ(f(1.0, 3.0), 7.0);
}

TEST(AstaraParser, VariableOrderMatters) {
  astara::Function f("x - y", {"x", "y"});
  EXPECT_DOUBLE_EQ(f(5.0, 2.0), 3.0);
}

TEST(AstaraParser, sin) {
  astara::Function f("sin(0)", {});
  EXPECT_DOUBLE_EQ(f(), 0.0);
}

TEST(AstaraParser, cos) {
  astara::Function f("cos(0)", {});
  EXPECT_DOUBLE_EQ(f(), 1.0);
}

TEST(AstaraParser, sqrt) {
  astara::Function f("sqrt(9)", {});
  EXPECT_DOUBLE_EQ(f(), 3.0);
}

TEST(AstaraParser, Functions) {
  astara::Function f("sqrt(abs(-9))", {});
  EXPECT_DOUBLE_EQ(f(), 3.0);
}

TEST(AstaraParser, FunctionWithVariable) {
  astara::Function f("exp(x)", {"x"});
  EXPECT_DOUBLE_EQ(f(0.0), 1.0);
}

TEST(AstaraParser, ComplexExpression) {
  astara::Function f("sin(x)/exp(-sin(x))+(y-x/cos(x))^2 - sqrt(x)",
                     {"x", "y"});
  double result = f(3.0, 10.0);
  double x = 3;
  double y = 10;
  EXPECT_EQ(result, sin(x) / exp(-sin(x)) +
                        (y - x / cos(x)) * (y - x / cos(x)) - sqrt(x));
}

TEST(AstaraParser, WhitespaceTolerance) {
  astara::Function f("  2 +   3 *   ( 4 + 1 ) ", {});
  EXPECT_DOUBLE_EQ(f(), 17.0);
}

TEST(AstaraParser, UnknownVariableThrows) {
  astara::Function f("x + 1", {});
  EXPECT_THROW(f(), std::runtime_error);
}

TEST(AstaraParser, FunctionThrows) {
  astara::Function f("foo(2)", {});
  EXPECT_THROW(f(), std::out_of_range);
}

TEST(AstaraParser, SyntaxErrorThrows) {
  astara::Function f("2 + * 3", {});
  EXPECT_THROW(f(), std::runtime_error);
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}