#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#include "Function.h"
#include "Reactor.h"

using namespace astara;

TEST(ReactorConstructorTest, ThrowsWhenGroupCountDoesNotMatchConstants) {
  EXPECT_THROW(Reactor reactor(3, {0.1, 0.2}, {0.1, 0.2}, 1e-9),
               std::runtime_error);
}

TEST(ReactorConstructorTest, DoesNotThrowWhenGroupCountMatchesConstants) {
  EXPECT_NO_THROW(Reactor reactor(2, {0.1, 0.2}, {0.1, 0.2}, 1e-9));
}

TEST(ReactorFunctionSetterTest, SettingAllFunctionsDoesNotThrow) {
  Reactor reactor(2, {0.1, 0.2}, {0.1, 0.2}, 1e-9);
  EXPECT_NO_THROW({
    reactor.setHeatTransferCoefficient(new Function("x", {"x"}));
    reactor.setFuelSpecificHeatFunction(new Function("1000 + T", {"T"}));
    reactor.setFuelTemperatureCoEfficientFunction(
        new Function("-0.01 * T", {"T"}));
    reactor.setModeratorTemperatureCoEfficientFunction(
        new Function("-0.002 * T", {"T"}));
  });
}
