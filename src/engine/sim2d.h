#pragma once
#include "mac/mac2d.h"
#include "solver/ivp.h"

class Simulation2D {
public:
  Simulation2D(int nx, int ny, float dx, float dy);
  ~Simulation2D() = default;

  void step(float dt);

private:
  MAC2D mac;
  MAC2D mac_next;

  int nx, ny;
  float dx, dy;

  IVPSolverType solver = IVPSolverType::EULER;

  void _advect_velocities(float dt);
  void _advect_u(float dt);
};
