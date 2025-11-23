#include "linear_ode_solver.h"

template <typename Q, typename... Params>
LinearODESolver<Q, Params...>::LinearODESolver(
    Q (*derivative)(Q, double, Params...), Q init_state, IVPSolverType solver)
    : _dq_dt(derivative), _state(init_state), _solver(solver) {}

template <typename Q, typename... Params>
void LinearODESolver<Q, Params...>::solve(double dt, Params... values) {
  _state = IVPSolvers<Q, Params...>::solveIVP(_solver, _state, _t, dt, _dq_dt,
                                              values...);
  _t += dt;
  _steps += 1;
}
