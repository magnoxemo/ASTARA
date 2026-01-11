#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#include "Reactor.h"

TEST(ReactorConstructorTest, ThrowsWhenGroupCountDoesNotMatchConstants) {
  unsigned int n_groups = 3;
  std::vector<double> delayed_constants = {0.1, 0.2};

  EXPECT_THROW(
      { astara::Reactor reactor(n_groups, delayed_constants); },
      std::runtime_error);
}

TEST(ReactorConstructorTest, DoesNotThrowWhenGroupCountMatchesConstants) {
  unsigned int n_groups = 2;
  std::vector<double> delayed_constants = {0.1, 0.2};

  EXPECT_NO_THROW({ astara::Reactor reactor(n_groups, delayed_constants); });
}
