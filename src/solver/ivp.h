#pragma once
#include <stdexcept>

enum IVPSolverType { EULER, MIDPOINT, RK4 };

template <typename Q, typename... Params> struct IVPSolvers {
  static Q eulerSolver(Q q0, float t0, float h, Q (*dq_dt)(Q, float, Params...),
                       Params... values) {
    return q0 + h * dq_dt(q0, t0, values...);
  }

  static Q midpointSolver(Q q0, float t0, float h,
                          Q (*dq_dt)(Q, float, Params...), Params... values) {
    Q k1 = dq_dt(q0, t0, values...);
    Q k2 = dq_dt(q0 + 0.5 * h * k1, t0 + 0.5 * h, values...);
    return q0 + h * k2;
  }

  static Q rk4Solver(Q q0, float t0, float h, Q (*dq_dt)(Q, float, Params...),
                     Params... values) {
    // from the wikipedia formulation:
    // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    Q k1 = dq_dt(q0, t0, values...);
    Q k2 = dq_dt(q0 + 0.5 * h * k1, t0 + 0.5 * h, values...);
    Q k3 = dq_dt(q0 + 0.5 * h * k2, t0 + 0.5 * h, values...);
    Q k4 = dq_dt(q0 + h * k3, t0 + h, values...);
    return q0 + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
  }

  static Q solveIVP(IVPSolverType solver, Q q0, float t0, float h,
                    Q (*dq_dt)(Q, float, Params...), Params... values) {
    switch (solver) {
    case EULER:
      return IVPSolvers<Q, Params...>::eulerSolver(q0, t0, h, dq_dt, values...);
    case MIDPOINT:
      return IVPSolvers<Q, Params...>::midpointSolver(q0, t0, h, dq_dt,
                                                      values...);
    case RK4:
      return IVPSolvers<Q, Params...>::rk4Solver(q0, t0, h, dq_dt, values...);
    default:
      throw std::invalid_argument("Invalid IVP Solver");
    }
  }
};
