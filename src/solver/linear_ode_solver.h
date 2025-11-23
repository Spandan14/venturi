#pragma once

#include "ivp.h"

template <typename Q, typename... Params> class LinearODESolver {
public:
  LinearODESolver(Q (*derivative)(Q, double, Params...), Q init_state,
                  IVPSolverType solver);
  void solve(double dt, Params... values);

  [[nodiscard]] Q getState() { return _state; }
  void setState(Q state) { _state = state; }

private:
  Q (*_dq_dt)(Q, double, Params...);

  int _steps = 0;
  double _t = 0.0;
  Q _state{};

  IVPSolverType _solver;
};
