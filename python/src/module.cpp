/**
 * @file   module.cpp
 * @brief  Top-level nanobind module entry point for `_astara`.
 *
 * Each `bind_<module>` function lives in its own translation unit, mirroring
 * the C++ library's own `src/<module>/` layout. `astara/__init__.py`
 * re-exports the symbols this registers under a clean `astara.*` namespace.
 */

#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_core(nb::module_& m);
void bind_props(nb::module_& m);
void bind_reactor(nb::module_& m);
void bind_pump(nb::module_& m);
void bind_pressurizer(nb::module_& m);
void bind_sg(nb::module_& m);
void bind_control(nb::module_& m);
void bind_primary(nb::module_& m);

NB_MODULE(_astara, m) {
    m.doc() = "Python bindings for the ASTARA PWR primary-loop simulator";
    bind_core(m);
    bind_props(m);
    bind_reactor(m);
    bind_pump(m);
    bind_pressurizer(m);
    bind_sg(m);
    bind_control(m);
    bind_primary(m);
}
