#include "astara/primary/PrimaryLoop.hpp"
#include "astara/control/ReactorController.hpp"
#include "astara/control/PressurizerController.hpp"
#include "astara/control/ThreeElementController.hpp"

#include <nanobind/nanobind.h>

#include <memory>
#include <utility>

namespace nb = nanobind;
using namespace astara::primary;

void bind_primary(nb::module_& m) {
    // PrimaryLoop's real C++ constructor takes plain std::unique_ptr<T>
    // (std::default_delete) component arguments, which nanobind cannot
    // transfer ownership into from a Python-constructed instance (those live
    // in memory embedded inside the Python object, not behind `new`; see
    // nanobind's ownership.rst on std::unique_ptr<T, nb::deleter<T>>). The
    // components are plain value types (no unique_ptr members of their own),
    // so instead of fighting that, the binding copies each Python-side
    // component into a freshly heap-allocated instance that PrimaryLoop can
    // safely own -- the same ownership-transfer semantics as the C++ API,
    // implemented via copy rather than a raw pointer handoff.
    nb::class_<PrimaryLoop>(m, "PrimaryLoop")
        .def(nb::new_([](astara::reactor::Reactor reactor,
                         astara::sg::AliSteamGenerator sg,
                         astara::pump::HomologousPump pump,
                         astara::pressurizer::Pressurizer pressurizer) {
                 return new PrimaryLoop(
                     std::make_unique<astara::reactor::Reactor>(std::move(reactor)),
                     std::make_unique<astara::sg::AliSteamGenerator>(std::move(sg)),
                     std::make_unique<astara::pump::HomologousPump>(std::move(pump)),
                     std::make_unique<astara::pressurizer::Pressurizer>(std::move(pressurizer)));
             }),
             nb::arg("reactor"), nb::arg("sg"), nb::arg("pump"), nb::arg("pressurizer"))
        .def("reactor",
             static_cast<astara::reactor::Reactor& (PrimaryLoop::*)()>(&PrimaryLoop::reactor),
             nb::rv_policy::reference_internal)
        .def("steam_generator",
             static_cast<astara::sg::AliSteamGenerator& (PrimaryLoop::*)()>(&PrimaryLoop::steamGenerator),
             nb::rv_policy::reference_internal)
        .def("pump",
             static_cast<astara::pump::HomologousPump& (PrimaryLoop::*)()>(&PrimaryLoop::pump),
             nb::rv_policy::reference_internal)
        .def("pressurizer",
             static_cast<astara::pressurizer::Pressurizer& (PrimaryLoop::*)()>(&PrimaryLoop::pressurizer),
             nb::rv_policy::reference_internal)
        // Same copy-into-a-fresh-unique_ptr rationale as the constructor above.
        .def("set_reactor_controller",
             [](PrimaryLoop& self, astara::control::ReactorController c) {
                 self.setReactorController(
                     std::make_unique<astara::control::ReactorController>(std::move(c)));
             })
        .def("set_pressurizer_controller",
             [](PrimaryLoop& self, astara::control::PressurizerController c) {
                 self.setPressurizerController(
                     std::make_unique<astara::control::PressurizerController>(std::move(c)));
             })
        .def("set_feedwater_controller",
             [](PrimaryLoop& self, astara::control::ThreeElementController c) {
                 self.setFeedwaterController(
                     std::make_unique<astara::control::ThreeElementController>(std::move(c)));
             })
        .def("has_reactor_controller",     &PrimaryLoop::hasReactorController)
        .def("has_pressurizer_controller", &PrimaryLoop::hasPressurizerController)
        .def("has_feedwater_controller",   &PrimaryLoop::hasFeedwaterController)
        .def("time_step",     &PrimaryLoop::timeStep)
        .def("time_seconds",  &PrimaryLoop::timeSeconds)
        .def("is_consistent", &PrimaryLoop::isConsistent, nb::arg("tol_K") = 5.0);
}
