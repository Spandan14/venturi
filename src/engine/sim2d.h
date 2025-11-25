#pragma once

#include "mac/mac2d.h"
#include "solver/ivp.h"
#include "solver/pressure_solvers.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

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
  void initialize_solids(const std::function<bool(int, int)> &initializer);

private:
  MAC2D mac;
  MAC2D mac_next;

  int nx, ny;
  float dx, dy;

  IVPSolverType ivp_solver = IVPSolverType::EULER;
  PressureSolverType pressure_solver = PressureSolverType::ICCG;

  void _apply_forces(float dt);

  void _advect_velocities(float dt);
  void _advect_u(float dt);
  void _advect_v(float dt);

  void _pressure_solve(float dt);
  void _project_velocities(const Eigen::VectorXd &pressure,
                           const std::vector<int> &fluid_idx, float dt);

  Eigen::VectorXd _build_divergence(const std::vector<int> &fluid_idx,
                                    int fluid_count);
  Eigen::SparseMatrix<double>
  _build_pressure_matrix(const std::vector<int> &fluid_idx, int fluid_count);

  void _advect_cell_data(float dt);
};
