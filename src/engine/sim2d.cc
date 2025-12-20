#include "sim2d.h"
#include "solver/iccg.h"
#include <fstream>
#include <iostream>
#include <utils/physical_consts.h>

Simulation2D::Simulation2D(int nx, int ny, float dx, float dy)
    : mac(nx, ny, dx, dy), mac_next(nx, ny, dx, dy), nx(nx), ny(ny), dx(dx),
      dy(dy) {}

void Simulation2D::step(float dt) {
  float total_density = 0.0f;
  for (float d : mac.densities) {
    total_density += d;
  }

  std::cout << "Total Density: " << total_density << std::endl;
  std::cout << "Max density: "
            << *std::max_element(mac.densities.begin(), mac.densities.end())
            << std::endl;

  _apply_forces(dt);      // will apply forces in mac
  _advect_velocities(dt); // will store advected velocities in mac_next

  _pressure_solve(dt);

  _advect_cell_data(dt);

  mac_next.current_time += dt;

  _apply_flows();
  _apply_flow_ratios();

  // total_density = 0.0f;
  // for (float d : mac_next.densities) {
  //   total_density += d;
  // }
  // std::cout << "Total Density After Step: " << total_density << std::endl;
  //
  // everything has now ended up in mac_next, so we swap
  std::swap(mac, mac_next);
}

void Simulation2D::initialize_vel_u(
    const std::function<float(int, int)> &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) {
      mac.u[mac.u_idx(i, j)] = initializer(i, j);
      mac_next.u[mac.u_idx(i, j)] = initializer(i, j);
    }
  }
}

void Simulation2D::initialize_vel_v(
    const std::function<float(int, int)> &initializer) {
  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.v[mac.v_idx(i, j)] = initializer(i, j);
      mac_next.v[mac.v_idx(i, j)] = initializer(i, j);
    }
  }
}

void Simulation2D::initialize_density(const DensityInitializer &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.densities[mac.idx(i, j)] = initializer(i, j, 0);
      mac_next.densities[mac.idx(i, j)] = initializer(i, j, 0);
    }
  }
}

void Simulation2D::initialize_forces(const ForceInitializer &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.forces[mac.idx(i, j)] = initializer(i, j, 0);
      mac_next.forces[mac.idx(i, j)] = initializer(i, j, 0);
    }
  }
}

void Simulation2D::initialize_solids(const SolidInitializer &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.is_solid[mac.idx(i, j)] = initializer(i, j, 0);
      mac_next.is_solid[mac.idx(i, j)] = initializer(i, j, 0);

      if (mac.is_solid[mac.idx(i, j)]) {
        // solid cell, set velocities to 0
        mac.u[mac.u_idx(i, j)] = 0.0f;
        mac.u[mac.u_idx(i + 1, j)] = 0.0f;
        mac.v[mac.v_idx(i, j)] = 0.0f;
        mac.v[mac.v_idx(i, j + 1)] = 0.0f;
        mac_next.u[mac.u_idx(i, j)] = 0.0f;
        mac_next.u[mac.u_idx(i + 1, j)] = 0.0f;
        mac_next.v[mac.v_idx(i, j)] = 0.0f;
        mac_next.v[mac.v_idx(i, j + 1)] = 0.0f;

        mac.densities[mac.idx(i, j)] = 0.00f;
        mac_next.densities[mac.idx(i, j)] = 0.00f;
      }
    }
  }
}

void Simulation2D::initialize_flows(const FlowGenerator &generator) {
  flow_generator = generator;
}

void Simulation2D::initialize_flow_ratios(const FlowRatioGenerator &generator) {
  flow_ratio_generator = generator;
}

void Simulation2D::_apply_forces(float dt) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) { // u-grid
      float force = 0.0f;
      if (i > 0) {
        force += mac.forces[mac.idx(i - 1, j)][0];
      }
      if (i < nx) {
        force += mac.forces[mac.idx(i, j)][0];
      }
      force *= 0.5f;

      mac.u[mac.u_idx(i, j)] += force * dt;
    }
  }

  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx; ++i) { // v-grid
      float force = 0.0f;
      if (j > 0) {
        force += mac.forces[mac.idx(i, j - 1)][1];
      }
      if (j < ny) {
        force += mac.forces[mac.idx(i, j)][1];
      }
      force *= 0.5f;

      mac.v[mac.v_idx(i, j)] += force * dt;
    }
  }
}

void Simulation2D::_apply_flows() {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      float added_density = flow_generator(i, j, mac_next.current_time);
      mac_next.densities[mac_next.idx(i, j)] += added_density;
      mac_next.densities[mac_next.idx(i, j)] =
          std::max(0.0f, mac_next.densities[mac_next.idx(i, j)]);
    }
  }
}

void Simulation2D::_apply_flow_ratios() {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      float ratio = flow_ratio_generator(i, j, mac_next.current_time);
      mac_next.densities[mac_next.idx(i, j)] *= ratio;
      mac_next.densities[mac_next.idx(i, j)] =
          std::max(0.0f, mac_next.densities[mac_next.idx(i, j)]);
    }
  }
}

// advects velocities using the semi-Lagrangian methods
// looks at every destination velocity grid point by axis, solves for the
// original, and advects new velocity
void Simulation2D::_advect_velocities(float dt) {
  _advect_u(dt);
  _advect_v(dt);
}

void Simulation2D::_advect_u(float dt) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) { // u-grid
      float x = i * mac.dx;
      float y = (j + 0.5f) * mac.dy;

      vec2d x_dest = vec2d(x, y);

      // solve for original position
      vec2d x_orig = IVPSolvers<vec2d, MAC2D &>::solveIVP(
          ivp_solver, x_dest, mac.current_time, dt, MAC2D::dx_vel_dt, mac);

      // interpolate the velocity at the original position
      mac_next.u[mac.u_idx(i, j)] = mac.vel_u(x_orig);
    }
  }
}

void Simulation2D::_advect_v(float dt) {
  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx; ++i) { // v-grid
      float x = (i + 0.5f) * mac.dx;
      float y = j * mac.dy;

      vec2d x_dest = vec2d(x, y);

      // solve for original position
      vec2d x_orig = IVPSolvers<vec2d, MAC2D &>::solveIVP(
          ivp_solver, x_dest, mac.current_time, dt, MAC2D::dx_vel_dt, mac);

      // interpolate the velocity at the original position
      mac_next.v[mac.v_idx(i, j)] = mac.vel_v(x_orig);
    }
  }
}

void Simulation2D::_pressure_solve(float dt) {
  std::vector<int> fluid_idx(nx * ny, -1);
  int fluid_count = 0;

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (mac.get_cell_type(i, j) == CellType::FLUID) {
        fluid_idx[mac.idx(i, j)] = fluid_count;
        fluid_count++;
      }
    }
  }

  // std::cout << "Fluid Cell Count: " << fluid_count << std::endl;

  if (fluid_count == 0) {
    // no fluid cells, nothing to do
    return;
  }

  Eigen::VectorXd divergence = _build_divergence(fluid_idx, fluid_count);
  Eigen::SparseMatrix<double> A =
      _build_pressure_matrix(fluid_idx, fluid_count);

  Eigen::VectorXd pressure(fluid_count);

  std::ofstream divergence_file("divergence_" +
                                std::to_string(mac.current_time) + ".vec");
  // divergence_file << divergence;
  for (int i = 0; i < divergence.size(); ++i) {
    divergence_file << divergence(i) << "\n";
  }
  divergence_file.close();

  std::ofstream pressure_file("pressure_matrix_" +
                              std::to_string(mac.current_time) + ".mtx");
  // pressure_file << A;
  // pressure_file.close();
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      pressure_file << it.row() << " " << it.col() << " " << it.value() << "\n";
    }
  }
  pressure_file.close();

  switch (pressure_solver) {

  case PressureSolverType::ICCG: {
    ICCGSolver iccg_solver = ICCGSolver(1000, 1e-6f);
    pressure = iccg_solver.solve(divergence, A);
    break;
  }
  default:
    throw std::runtime_error("[SIMULATION] Unsupported pressure solver type.");
  }

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int c_idx = mac.idx(i, j);
      if (fluid_idx[c_idx] != -1) {
        mac_next.pressures[c_idx] = pressure[fluid_idx[c_idx]];
      } else {
        mac_next.pressures[c_idx] = 0.0f;
      }
    }
  }

  _project_velocities(pressure, fluid_idx, dt);
}

void Simulation2D::_project_velocities(const Eigen::VectorXd &pressure,
                                       const std::vector<int> &fluid_idx,
                                       float dt) {
  // u faces
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) { // internal faces only
      if (i == 0 || i == nx) {
        mac_next.u[mac.u_idx(i, j)] = 0.0f;
        continue;
      }

      int c_idx_left = mac.idx(i - 1, j);
      int c_idx_right = mac.idx(i, j);

      bool left_valid = (fluid_idx[c_idx_left] != -1);
      bool right_valid = (fluid_idx[c_idx_right] != -1);

      bool left_non_solid = mac_next.get_cell_type(i - 1, j) != CellType::SOLID;
      bool right_non_solid = mac_next.get_cell_type(i, j) != CellType::SOLID;

      // if (left_valid || right_valid) {
      if (left_non_solid && right_non_solid) {
        if (!left_valid && !right_valid) {
          mac_next.u[mac.u_idx(i, j)] = 0.0f;
          continue; // both are air
        }

        // air always has 0 pressure
        float p_left = left_valid ? pressure[fluid_idx[c_idx_left]] : 0.0f;
        float p_right = right_valid ? pressure[fluid_idx[c_idx_right]] : 0.0f;

        // average density
        float rho;
        if (left_valid && right_valid) {
          rho = 0.5f * (mac.densities[c_idx_left] + mac.densities[c_idx_right]);
        } else if (left_valid) {
          rho = mac.densities[c_idx_left];
        } else {
          rho = mac.densities[c_idx_right];
        }

        float correction = (p_right - p_left) * (dt / rho);
        correction /= mac.dx;

        mac_next.u[mac.u_idx(i, j)] -= correction;
      } else {
        // one side is solid, set to 0
        mac_next.u[mac.u_idx(i, j)] = 0.0f;
      }
      //   mac_next.u[mac.u_idx(i, j)] = 0.0f;
      // }
    }
  }

  // v faces
  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (j == 0 || j == ny) {
        mac_next.v[mac.v_idx(i, j)] = 0.0f;
        continue;
      }

      int c_idx_bottom = mac.idx(i, j - 1);
      int c_idx_top = mac.idx(i, j);

      bool bottom_valid = (fluid_idx[c_idx_bottom] != -1);
      bool top_valid = (fluid_idx[c_idx_top] != -1);

      bool bottom_non_solid =
          mac_next.get_cell_type(i, j - 1) != CellType::SOLID;
      bool top_non_solid = mac_next.get_cell_type(i, j) != CellType::SOLID;

      if (top_non_solid && bottom_non_solid) {
        if (!top_valid && !bottom_valid) {
          mac_next.v[mac.v_idx(i, j)] = 0.0f;
          continue; // both are air
        }

        // air always has 0 pressure
        float p_top = top_valid ? pressure[fluid_idx[c_idx_top]] : 0.0f;
        float p_bottom =
            bottom_valid ? pressure[fluid_idx[c_idx_bottom]] : 0.0f;

        // average density
        float rho;
        if (top_valid && bottom_valid) {
          rho = 0.5f * (mac.densities[c_idx_top] + mac.densities[c_idx_bottom]);
        } else if (top_valid) {
          rho = mac.densities[c_idx_top];
        } else {
          rho = mac.densities[c_idx_bottom];
        }

        float correction = (p_top - p_bottom) * (dt / rho);
        correction /= mac.dy;

        // std::cout << "Original v at (" << i << ", " << j
        //           << "): " << mac_next.v[mac.v_idx(i, j)] << std::endl;
        // std::cout << "p_top: " << p_top << ", p_bottom: " << p_bottom
        //           << ", rho: " << rho << ", dt: " << dt << std::endl;
        // std::cout << "Correction: " << correction << std::endl;

        mac_next.v[mac.v_idx(i, j)] -= correction;
      } else {
        // one side is solid, set to 0
        mac_next.v[mac.v_idx(i, j)] = 0.0f;
      }
      //   mac_next.v[mac.v_idx(i, j)] = 0.0f;
      // }
    }
  }
}

Eigen::VectorXd
Simulation2D::_build_divergence(const std::vector<int> &fluid_idx,
                                int fluid_count) {
  Eigen::VectorXd divergence(fluid_count);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int c_idx = mac.idx(i, j);
      if (fluid_idx[c_idx] == -1) {
        continue; // not a fluid cell
      }

      float div = 0.0f;

      // NOTE: KEY to use new vel here
      // u velocities
      // div += mac_next.u[mac.u_idx(i + 1, j)]; // right face
      div += mac.get_cell_type(i + 1, j) == CellType::SOLID
                 ? 0.0f
                 : mac_next.u[mac.u_idx(i + 1, j)]; // right face

      div -= mac.get_cell_type(i - 1, j) == CellType::SOLID
                 ? 0.0f
                 : mac_next.u[mac.u_idx(i, j)]; // left face

      // v velocities
      div += mac.get_cell_type(i, j + 1) == CellType::SOLID
                 ? 0.0f
                 : mac_next.v[mac.v_idx(i, j + 1)]; // top face

      div -= mac.get_cell_type(i, j - 1) == CellType::SOLID
                 ? 0.0f
                 : mac_next.v[mac.v_idx(i, j)]; // bottom face

      divergence[fluid_idx[c_idx]] = -div / mac.dx; // negative for RHS
    }
  }

  return divergence;
}

Eigen::SparseMatrix<double>
Simulation2D::_build_pressure_matrix(const std::vector<int> &fluid_idx,
                                     int fluid_count) {
  std::vector<Eigen::Triplet<double>> triplet_list;

  float inv_dx2 = 1.0f / (mac.dx * mac.dx);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int c_idx = mac.idx(i, j);
      if (fluid_idx[c_idx] == -1) {
        continue; // not a fluid cell
      }

      int mat_row = fluid_idx[c_idx];
      int neighbor_count = 0;

      // check neighbors
      const std::vector<std::pair<int, int>> neighbors = {
          {1, 0}, {-1, 0}, {0, 1}, {0, -1}};

      for (const auto &[di, dj] : neighbors) {
        int ni = i + di;
        int nj = j + dj;

        if (ni >= 0 && ni < nx && nj >= 0 && nj < ny) {
          int n_idx = mac.idx(ni, nj);
          if (fluid_idx[n_idx] != -1) {
            // fluid neighbor
            triplet_list.emplace_back(mat_row, fluid_idx[n_idx], -inv_dx2);
          }

          if (mac.get_cell_type(ni, nj) != CellType::SOLID) {
            // non-solid neighbor
            neighbor_count++;
          }
        }
      }

      // diagonal entry
      triplet_list.emplace_back(mat_row, mat_row,
                                static_cast<double>(neighbor_count) * inv_dx2);
    }
  }

  Eigen::SparseMatrix<double> A(fluid_count, fluid_count);
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());
  return A;
}

void Simulation2D::_advect_cell_data(float dt) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) { // cell-centered
      if (mac.get_cell_type(i, j) == CellType::SOLID) {
        // solid cell, skip
        mac_next.densities[mac.idx(i, j)] = mac.densities[mac.idx(i, j)];
        continue;
      }

      float x = (i + 0.5f) * mac.dx;
      float y = (j + 0.5f) * mac.dy;

      vec2d x_dest = vec2d(x, y);

      // solve for original position
      vec2d x_orig = IVPSolvers<vec2d, MAC2D &>::solveIVP(
          ivp_solver, x_dest, mac.current_time, dt, MAC2D::dx_vel_dt, mac_next);

      // interpolate the cell data at the original position
      // for now, just density

      if (mac.is_position_solid(x_orig)) {
        // std::cout << "Warning: Density advected into solid at cell (" << i
        //           << ", " << j << ")\n";
        //           << std::endl;
        x_orig = mac.nonsolid_projection(x_orig, x_dest,
                                         nonsolid_projection_iterations);
        // std::cout << "assigned density: " << mac.density(x_orig) <<
        // std::endl;
      }

      mac_next.densities[mac.idx(i, j)] = mac.density(x_orig);
    }
  }
}
