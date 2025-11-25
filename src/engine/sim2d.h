#pragma once
#include "mac/mac2d.h"
#include "solver/ivp.h"

class Simulation2D {
public:
  Simulation2D(int nx, int ny, float dx, float dy);
  ~Simulation2D() = default;

  void step(float dt);
  [[nodiscard]] const MAC2D &get_mac() const { return mac; }

  void initialize_vel_u(const std::function<float(int, int)> &initializer);
  void initialize_vel_v(const std::function<float(int, int)> &initializer);
  void initialize_density(const std::function<float(int, int)> &initializer);
  void initialize_forces(const std::function<vec2d(int, int)> &initializer);

private:
  MAC2D mac;
  MAC2D mac_next;

  int nx, ny;
  float dx, dy;

  IVPSolverType solver = IVPSolverType::EULER;

  void _advect_velocities(float dt);
  void _advect_u(float dt);
  void _advect_v(float dt);

  void _advect_cell_data(float dt);

  void _apply_forces(float dt);
};
