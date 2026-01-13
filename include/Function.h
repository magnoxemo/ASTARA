#ifndef ASTARA_PARSER_H
#define ASTARA_PARSER_H

#include <cctype>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace astara {

class Function {
public:
  Function(const std::string &expr,
           const std::vector<std::string> &variable_names);

  // Call operator: x=values[0], y=values[1], etc.
  template <typename... Args> double operator()(Args... args) {
    _values = {static_cast<double>(args)...};

    if (_values.size() != _variables.size()) {
      throw std::runtime_error("Incorrect number of arguments");
    }

    _pos = 0;
    double result = parseExpression();
    skipWhitespace();

    if (_pos != _expression.size()) {
      throw std::runtime_error("Unexpected trailing characters");
    }

    return result;
  }

private:
  // Grammar:
  // expression = term { (+|-) term }
  // term       = factor { (*|/) factor }
  // factor     = unary [ ^ factor ]
  // unary      = (+|-) unary | primary
  // primary    = number | variable | function | '(' expression ')'

  double parseExpression();
  double parseTerm();
  double parseFactor();
  double parseUnary();
  double parsePrimary();
  double parseNumber();
  std::string parseIdentifier();
  double getVariable(const std::string &name);
  double callFunction(const std::string &name, double arg);
  bool match(char c);
  char peek() const;
  void skipWhitespace();

  // variables
  const std::string _expression;
  size_t _pos;
  const std::vector<std::string> _variables;
  std::vector<double> _values;

  std::unordered_map<std::string, std::function<double(double)>> _functions = {
      {"sin", static_cast<double (*)(double)>(std::sin)},
      {"cos", static_cast<double (*)(double)>(std::cos)},
      {"tan", static_cast<double (*)(double)>(std::tan)},
      {"log", static_cast<double (*)(double)>(std::log)},
      {"sqrt", static_cast<double (*)(double)>(std::sqrt)},
      {"exp", static_cast<double (*)(double)>(std::exp)},
      {"abs", static_cast<double (*)(double)>(std::fabs)}};
};

} // namespace astara

#endif
