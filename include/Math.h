#ifndef ASTARA_MATH_H
#define ASTARA_MATH_H

// main math functions and solvers for the framework
namespace astara{

    template <typename T>
    T predictor_corrector_integrator(T (*differential_func)(), T time_step);

    template <typename T>
    T predictor_integrator(T (*differential_func)(), T time_step);
}

#endif //ASTARA_MATH_H


