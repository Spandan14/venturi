#pragma once

#include "engine/sim.h"
#include "mac/mac3d.h"
#include "solver/ivp.h"
#include "solver/pressure_solvers.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

class Simulation3D : public Simulation<3> {
public:
  Simulation3D(int nx, int ny, int nz, float dx, float dy, float dz);
  ~Simulation3D() = default;

  void step(float dt) override;
  [[nodiscard]] const MAC3D &get_mac() const { return mac; }

  void initialize_density(const DensityInitializer &initializer) override;
  void initialize_forces(const ForceInitializer &initializer) override;
  void initialize_solids(const SolidInitializer &initializer) override;

  void initialize_flows(const FlowGenerator &generator) override;
  void initialize_flow_ratios(const FlowRatioGenerator &generator) override;

private:
  MAC3D mac;
  MAC3D mac_next;

  int nx, ny, nz;
  float dx, dy, dz;

  IVPSolverType ivp_solver = IVPSolverType::RK4;
  PressureSolverType pressure_solver = PressureSolverType::ICCG;
  int nonsolid_projection_iterations = 20;

  void _apply_forces(float dt);

  void _advect_velocities(float dt);
  void _advect_u(float dt);
  void _advect_v(float dt);
  void _advect_w(float dt);

  void _pressure_solve(float dt);
  void _project_velocities(const Eigen::VectorXd &pressure,
                           const std::vector<int> &fluid_idx, float dt);

  Eigen::VectorXd _build_divergence(const std::vector<int> &fluid_idx,
                                    int fluid_count);
  Eigen::SparseMatrix<double>
  _build_pressure_matrix(const std::vector<int> &fluid_idx, int fluid_count);

  void _advect_cell_data(float dt);

  FlowGenerator flow_generator;
  void _apply_flows();

  FlowRatioGenerator flow_ratio_generator;
  void _apply_flow_ratios();

  SolidInitializer solid_initializer;
  void _apply_solids(const SolidInitializer &initializer);
};
