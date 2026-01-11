#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#include "Parser.h"
#include "Reactor.h"

using namespace astara;


TEST(ReactorConstructorTest, ThrowsWhenGroupCountDoesNotMatchConstants) {
  EXPECT_THROW(Reactor reactor(3, {0.1, 0.2}), std::runtime_error);
}

TEST(ReactorConstructorTest, DoesNotThrowWhenGroupCountMatchesConstants) {
  EXPECT_NO_THROW(Reactor reactor(2, {0.1, 0.2}));
}

// ------------------------
// Heat transfer coefficient setter
// ------------------------
TEST(ReactorFunctionSetterTest, HeatTransferCoefficientThrowsOnNullptr) {
  Reactor reactor(1, {0.1});
  EXPECT_THROW(reactor.setHeatTransferCoefficient(nullptr), std::runtime_error);
}

TEST(ReactorFunctionSetterTest, HeatTransferCoefficientAcceptsValidFunction) {
  Reactor reactor(1, {0.1});
  EXPECT_NO_THROW(
      reactor.setHeatTransferCoefficient(new Function("2.0 * x", {"x"})));
}

// ------------------------
// Fuel specific heat setter
// ------------------------
TEST(ReactorFunctionSetterTest, FuelSpecificHeatThrowsOnNullptr) {
  Reactor reactor(1, {0.1});
  EXPECT_THROW(reactor.setFuelSpecificHeatFunction(nullptr),
               std::runtime_error);
}

TEST(ReactorFunctionSetterTest, FuelSpecificHeatAcceptsValidFunction) {
  Reactor reactor(1, {0.1});
  EXPECT_NO_THROW(
      reactor.setFuelSpecificHeatFunction(new Function("1000 + T", {"T"})));
}

// ------------------------
// Fuel temperature coefficient setter
// ------------------------
TEST(ReactorFunctionSetterTest, FuelTemperatureCoefficientThrowsOnNullptr) {
  Reactor reactor(1, {0.1});
  EXPECT_THROW(reactor.setFuelTemperatureCoEfficientFunction(nullptr),
               std::runtime_error);
}

TEST(ReactorFunctionSetterTest,
     FuelTemperatureCoefficientAcceptsValidFunction) {
  Reactor reactor(1, {0.1});
  EXPECT_NO_THROW(reactor.setFuelTemperatureCoEfficientFunction(
      new Function("-0.01 * T", {"T"})));
}

// ------------------------
// Moderator temperature coefficient setter
// ------------------------
TEST(ReactorFunctionSetterTest,
     ModeratorTemperatureCoefficientThrowsOnNullptr) {
  Reactor reactor(1, {0.1});
  EXPECT_THROW(reactor.setModeratorTemperatureCoEfficientFunction(nullptr),
               std::runtime_error);
}

TEST(ReactorFunctionSetterTest,
     ModeratorTemperatureCoefficientAcceptsValidFunction) {
  Reactor reactor(1, {0.1});
  EXPECT_NO_THROW(reactor.setModeratorTemperatureCoEfficientFunction(
      new Function("-0.002 * T", {"T"})));
}

// ------------------------
// Replacement semantics
// ------------------------
TEST(ReactorFunctionSetterTest, ReplacingHeatTransferCoefficientDoesNotThrow) {
  Reactor reactor(1, {0.1});
  EXPECT_NO_THROW({
    reactor.setHeatTransferCoefficient(new Function("x", {"x"}));
    reactor.setHeatTransferCoefficient(new Function("2 * x", {"x"}));
  });
}

// ------------------------
// Multiple setters together
// ------------------------
TEST(ReactorFunctionSetterTest, SettingAllFunctionsDoesNotThrow) {
  Reactor reactor(1, {0.1});
  EXPECT_NO_THROW({
    reactor.setHeatTransferCoefficient(new Function("x", {"x"}));
    reactor.setFuelSpecificHeatFunction(new Function("1000 + T", {"T"}));
    reactor.setFuelTemperatureCoEfficientFunction(
        new Function("-0.01 * T", {"T"}));
    reactor.setModeratorTemperatureCoEfficientFunction(
        new Function("-0.002 * T", {"T"}));
  });
}
