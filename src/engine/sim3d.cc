#include "sim3d.h"
#include "solver/iccg.h"
#include <iostream>
#include <utils/physical_consts.h>

Simulation3D::Simulation3D(int nx, int ny, int nz, float dx, float dy, float dz)
    : mac(nx, ny, nz, dx, dy, dz), mac_next(nx, ny, nz, dx, dy, dz), nx(nx),
      ny(ny), nz(nz), dx(dx), dy(dy), dz(dz) {}

void Simulation3D::step(float dt) {
  // float total_density = 0.0f;
  // for (float d : mac.densities) {
  //   total_density += d;
  // }
  //
  // std::cout << "Total Density: " << total_density << std::endl;
  // std::cout << "Max density: "
  //           << *std::max_element(mac.densities.begin(), mac.densities.end())
  //           << std::endl;
  //
  // _apply_solids(solid_initializer);
  _apply_forces(dt);
  _advect_velocities(dt);

  _pressure_solve(dt);

  _advect_cell_data(dt);

  mac_next.current_time += dt;

  _apply_flows();
  _apply_flow_ratios();

  // std::cout << "Simulation time: " << mac_next.current_time << " s"
  //           << std::endl;

  std::swap(mac, mac_next);
}

void Simulation3D::initialize_density(const DensityInitializer &initializer) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        mac.densities[mac.idx(i, j, k)] =
            initializer(i, j, k, mac.current_time);
        mac_next.densities[mac.idx(i, j, k)] =
            initializer(i, j, k, mac.current_time);
      }
    }
  }
}

void Simulation3D::initialize_forces(const ForceInitializer &initializer) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        mac.forces[mac.idx(i, j, k)] = initializer(i, j, k, mac.current_time);
        mac_next.forces[mac.idx(i, j, k)] =
            initializer(i, j, k, mac.current_time);
      }
    }
  }
}

void Simulation3D::initialize_solids(const SolidInitializer &initializer) {
  solid_initializer = initializer;

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        mac.is_solid[mac.idx(i, j, k)] = initializer(i, j, k, mac.current_time);
        mac_next.is_solid[mac.idx(i, j, k)] =
            initializer(i, j, k, mac.current_time);

        if (mac.is_solid[mac.idx(i, j, k)]) {
          // solid cell, set velocities to 0
          mac.u[mac.u_idx(i, j, k)] = 0.0f;
          mac.u[mac.u_idx(i + 1, j, k)] = 0.0f;
          mac.v[mac.v_idx(i, j, k)] = 0.0f;
          mac.v[mac.v_idx(i, j + 1, k)] = 0.0f;
          mac.w[mac.w_idx(i, j, k)] = 0.0f;
          mac.w[mac.w_idx(i, j, k + 1)] = 0.0f;

          mac_next.u[mac.u_idx(i, j, k)] = 0.0f;
          mac_next.u[mac.u_idx(i + 1, j, k)] = 0.0f;
          mac_next.v[mac.v_idx(i, j, k)] = 0.0f;
          mac_next.v[mac.v_idx(i, j + 1, k)] = 0.0f;
          mac_next.w[mac.w_idx(i, j, k)] = 0.0f;
          mac_next.w[mac.w_idx(i, j, k + 1)] = 0.0f;

          mac.densities[mac.idx(i, j, k)] = 0.00f;
          mac_next.densities[mac.idx(i, j, k)] = 0.00f;
        }
      }
    }
  }
}

void Simulation3D::initialize_flows(const FlowGenerator &generator) {
  flow_generator = generator;
}

void Simulation3D::initialize_flow_ratios(const FlowRatioGenerator &generator) {
  flow_ratio_generator = generator;
}

void Simulation3D::_apply_forces(float dt) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx + 1; ++i) { // u-grid
        float force = 0.0f;
        if (i > 0) {
          force += mac.forces[mac.idx(i - 1, j, k)][0];
        }
        if (i < nx) {
          force += mac.forces[mac.idx(i, j, k)][0];
        }
        force *= 0.5f;

        mac.u[mac.u_idx(i, j, k)] += force * dt;
      }
    }
  }

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny + 1; ++j) {
      for (int i = 0; i < nx; ++i) { // v-grid
        float force = 0.0f;
        if (j > 0) {
          force += mac.forces[mac.idx(i, j - 1, k)][1];
        }
        if (j < ny) {
          force += mac.forces[mac.idx(i, j, k)][1];
        }
        force *= 0.5f;

        mac.v[mac.v_idx(i, j, k)] += force * dt;
      }
    }
  }

  for (int k = 0; k < nz + 1; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) { // w-grid
        float force = 0.0f;
        if (k > 0) {
          force += mac.forces[mac.idx(i, j, k - 1)][2];
        }
        if (k < nz) {
          force += mac.forces[mac.idx(i, j, k)][2];
        }
        force *= 0.5f;

        mac.w[mac.w_idx(i, j, k)] += force * dt;
      }
    }
  }
}

void Simulation3D::_apply_flows() {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        float amount = flow_generator(i, j, k, mac_next.current_time);
        mac_next.densities[mac_next.idx(i, j, k)] += amount;
        mac_next.densities[mac_next.idx(i, j, k)] =
            std::max(0.0f, mac_next.densities[mac_next.idx(i, j, k)]);
      }
    }
  }
}

void Simulation3D::_apply_flow_ratios() {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        float ratio = flow_ratio_generator(i, j, k, mac_next.current_time);
        mac_next.densities[mac_next.idx(i, j, k)] *= ratio;
        mac_next.densities[mac_next.idx(i, j, k)] =
            std::max(0.0f, mac_next.densities[mac_next.idx(i, j, k)]);
      }
    }
  }
}

void Simulation3D::_apply_solids(const SolidInitializer &initializer) {
  std::vector<std::tuple<int, int, int>> previously_solid;
  float total_released_density = 0.0f;

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        bool was_solid = mac.is_solid[mac.idx(i, j, k)];

        mac.is_solid[mac.idx(i, j, k)] = initializer(i, j, k, mac.current_time);
        // mac_next.is_solid[mac.idx(i, j, k)] =
        //     initializer(i, j, k, mac.current_time);

        if (mac.is_solid[mac.idx(i, j, k)]) {
          // solid cell, set velocities to 0
          mac.u[mac.u_idx(i, j, k)] = 0.0f;
          mac.u[mac.u_idx(i + 1, j, k)] = 0.0f;
          mac.v[mac.v_idx(i, j, k)] = 0.0f;
          mac.v[mac.v_idx(i, j + 1, k)] = 0.0f;
          mac.w[mac.w_idx(i, j, k)] = 0.0f;
          mac.w[mac.w_idx(i, j, k + 1)] = 0.0f;

          float old_density = mac.densities[mac.idx(i, j, k)];
          total_released_density += old_density;

          mac.densities[mac.idx(i, j, k)] = 0.00f;

          // std::vector<std::tuple<float, float, float>> diffs = {
          //     {1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
          //     {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
          //
          // for (const auto &[di, dj, dk] : diffs) {
          //   int ni = i + di;
          //   int nj = j + dj;
          //   int nk = k + dk;
          //   if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 &&
          //       nk < nz) {
          //     if (!mac.is_solid[mac.idx(ni, nj, nk)]) {
          //       mac.densities[mac.idx(ni, nj, nk)] += old_density / 6.0f;
          //     }
          //   }
          // }

          // mac_next.u[mac.u_idx(i, j, k)] = 0.0f;
          // mac_next.u[mac.u_idx(i + 1, j, k)] = 0.0f;
          // mac_next.v[mac.v_idx(i, j, k)] = 0.0f;
          // mac_next.v[mac.v_idx(i, j + 1, k)] = 0.0f;
          // mac_next.w[mac.w_idx(i, j, k)] = 0.0f;
          // mac_next.w[mac.w_idx(i, j, k + 1)] = 0.0f;
          //
          // mac_next.densities[mac.idx(i, j, k)] = 0.00f;
        }

        if (was_solid && !mac.is_solid[mac.idx(i, j, k)]) {
          previously_solid.emplace_back(i, j, k);
        }
      }
    }
  }

  // redistribute density from previously solid cells
  for (const auto &[i, j, k] : previously_solid) {
    // float distribute_density =
    //     total_released_density / static_cast<float>(previously_solid.size());
    // mac.densities[mac.idx(i, j, k)] +=
    //     std::max(distribute_density, MIN_FLUID_DENSITY);
    // grab density from average of neighbors

    std::vector<std::tuple<float, float, float>> diffs = {
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};

    float accumulated_density = 0.0f;
    int valid_neighbors = 0;
    for (const auto &[di, dj, dk] : diffs) {
      int ni = i + di;
      int nj = j + dj;
      int nk = k + dk;
      if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
        if (!mac.is_solid[mac.idx(ni, nj, nk)]) {
          accumulated_density += mac.densities[mac.idx(ni, nj, nk)];
          valid_neighbors++;
        }
      }
    }

    if (valid_neighbors > 0) {
      mac.densities[mac.idx(i, j, k)] +=
          std::max(accumulated_density / static_cast<float>(valid_neighbors),
                   MIN_FLUID_DENSITY);
    } else {
      mac.densities[mac.idx(i, j, k)] += MIN_FLUID_DENSITY;
    }
  }
}

void Simulation3D::_advect_velocities(float dt) {
  _advect_u(dt);
  _advect_v(dt);
  _advect_w(dt);
}

void Simulation3D::_advect_u(float dt) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx + 1; ++i) { // u-grid
        float x = i * mac.dx;
        float y = (j + 0.5f) * mac.dy;
        float z = (k + 0.5f) * mac.dz;

        vec3d x_dest = vec3d(x, y, z);

        // solve for original position
        vec3d x_orig = IVPSolvers<vec3d, MAC3D &>::solveIVP(
            ivp_solver, x_dest, mac.current_time, dt, MAC3D::dx_vel_dt, mac);

        // interpolate the velocity at the original position
        mac_next.u[mac.u_idx(i, j, k)] = mac.vel_u(x_orig);
      }
    }
  }
}

void Simulation3D::_advect_v(float dt) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny + 1; ++j) {
      for (int i = 0; i < nx; ++i) { // v-grid
        float x = (i + 0.5f) * mac.dx;
        float y = j * mac.dy;
        float z = (k + 0.5f) * mac.dz;

        vec3d x_dest = vec3d(x, y, z);

        // solve for original position
        vec3d x_orig = IVPSolvers<vec3d, MAC3D &>::solveIVP(
            ivp_solver, x_dest, mac.current_time, dt, MAC3D::dx_vel_dt, mac);

        // interpolate the velocity at the original position
        mac_next.v[mac.v_idx(i, j, k)] = mac.vel_v(x_orig);
      }
    }
  }
}

void Simulation3D::_advect_w(float dt) {
  for (int k = 0; k < nz + 1; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) { // w-grid
        float x = (i + 0.5f) * mac.dx;
        float y = (j + 0.5f) * mac.dy;
        float z = k * mac.dz;

        vec3d x_dest = vec3d(x, y, z);

        // solve for original position
        vec3d x_orig = IVPSolvers<vec3d, MAC3D &>::solveIVP(
            ivp_solver, x_dest, mac.current_time, dt, MAC3D::dx_vel_dt, mac);

        // interpolate the velocity at the original position
        mac_next.w[mac.w_idx(i, j, k)] = mac.vel_w(x_orig);
      }
    }
  }
}

void Simulation3D::_pressure_solve(float dt) {
  std::vector<int> fluid_idx(nx * ny * nz, -1);
  int fluid_count = 0;

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (mac.get_cell_type(i, j, k) == CellType::FLUID) {
          fluid_idx[mac.idx(i, j, k)] = fluid_count;
          fluid_count++;
        }
      }
    }
  }

  if (fluid_count == 0) {
    return;
  }

  Eigen::VectorXd divergence = _build_divergence(fluid_idx, fluid_count);
  Eigen::SparseMatrix<double> A =
      _build_pressure_matrix(fluid_idx, fluid_count);

  Eigen::VectorXd pressure(fluid_count);

  switch (pressure_solver) {
  case PressureSolverType::ICCG: {
    ICCGSolver iccg_solver = ICCGSolver(1000, 1e-6f);
    pressure = iccg_solver.solve(divergence, A);
    break;
  }
  default:
    throw std::runtime_error("[SIMULATION] Unsupported pressure solver type.");
  }

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int c_idx = mac.idx(i, j, k);
        if (fluid_idx[c_idx] != -1) {
          mac_next.pressures[c_idx] = pressure[fluid_idx[c_idx]];
        } else {
          mac_next.pressures[c_idx] = 0.0f;
        }
      }
    }
  }

  _project_velocities(pressure, fluid_idx, dt);
}

void Simulation3D::_project_velocities(const Eigen::VectorXd &pressure,
                                       const std::vector<int> &fluid_idx,
                                       float dt) {
  // u faces
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx + 1; ++i) {
        if (i == 0 || i == nx) {
          mac_next.u[mac.u_idx(i, j, k)] = 0.0f;
          continue;
        }

        int c_idx_left = mac.idx(i - 1, j, k);
        int c_idx_right = mac.idx(i, j, k);

        bool left_valid = (fluid_idx[c_idx_left] != -1);
        bool right_valid = (fluid_idx[c_idx_right] != -1);

        bool left_non_solid =
            mac_next.get_cell_type(i - 1, j, k) != CellType::SOLID;
        bool right_non_solid =
            mac_next.get_cell_type(i, j, k) != CellType::SOLID;

        if (left_non_solid && right_non_solid) {
          if (!left_valid && !right_valid) {
            mac_next.u[mac.u_idx(i, j, k)] = 0.0f;
            continue; // both are air
          }

          float p_left = left_valid ? pressure[fluid_idx[c_idx_left]] : 0.0f;
          float p_right = right_valid ? pressure[fluid_idx[c_idx_right]] : 0.0f;

          float rho;
          if (left_valid && right_valid) {
            rho =
                0.5f * (mac.densities[c_idx_left] + mac.densities[c_idx_right]);
          } else if (left_valid) {
            rho = mac.densities[c_idx_left];
          } else {
            rho = mac.densities[c_idx_right];
          }

          float correction = (p_right - p_left) * (dt / rho);
          correction /= mac.dx;

          mac_next.u[mac.u_idx(i, j, k)] -= correction;
        } else {
          mac_next.u[mac.u_idx(i, j, k)] = 0.0f;
        }
      }
    }
  }

  // v faces
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny + 1; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (j == 0 || j == ny) {
          mac_next.v[mac.v_idx(i, j, k)] = 0.0;
          continue;
        }

        int c_idx_bottom = mac.idx(i, j - 1, k);
        int c_idx_top = mac.idx(i, j, k);

        bool bottom_valid = (fluid_idx[c_idx_bottom] != -1);
        bool top_valid = (fluid_idx[c_idx_top] != -1);

        bool bottom_non_solid =
            mac_next.get_cell_type(i, j - 1, k) != CellType::SOLID;
        bool top_non_solid = mac_next.get_cell_type(i, j, k) != CellType::SOLID;

        if (top_non_solid && bottom_non_solid) {
          if (!top_valid && !bottom_valid) {
            mac_next.v[mac.v_idx(i, j, k)] = 0.0f;
            continue; // both are air
          }

          float p_top = top_valid ? pressure[fluid_idx[c_idx_top]] : 0.0f;
          float p_bottom =
              bottom_valid ? pressure[fluid_idx[c_idx_bottom]] : 0.0f;

          float rho;
          if (top_valid && bottom_valid) {
            rho =
                0.5f * (mac.densities[c_idx_top] + mac.densities[c_idx_bottom]);
          } else if (top_valid) {
            rho = mac.densities[c_idx_top];
          } else {
            rho = mac.densities[c_idx_bottom];
          }

          float correction = (p_top - p_bottom) * (dt / rho);
          correction /= mac.dy;

          mac_next.v[mac.v_idx(i, j, k)] -= correction;
        } else {
          mac_next.v[mac.v_idx(i, j, k)] = 0.0f;
        }
      }
    }
  }

  // w faces
  for (int k = 0; k < nz + 1; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (k == 0 || k == nz) {
          mac_next.w[mac.w_idx(i, j, k)] = 0.0f;
          continue;
        }

        int c_idx_back = mac.idx(i, j, k - 1);
        int c_idx_front = mac.idx(i, j, k);

        bool back_valid = (fluid_idx[c_idx_back] != -1);
        bool front_valid = (fluid_idx[c_idx_front] != -1);

        bool back_non_solid =
            mac_next.get_cell_type(i, j, k - 1) != CellType::SOLID;
        bool front_non_solid =
            mac_next.get_cell_type(i, j, k) != CellType::SOLID;

        if (front_non_solid && back_non_solid) {
          if (!front_valid && !back_valid) {
            mac_next.w[mac.w_idx(i, j, k)] = 0.0f;
            continue; // both are air
          }

          float p_front = front_valid ? pressure[fluid_idx[c_idx_front]] : 0.0f;
          float p_back = back_valid ? pressure[fluid_idx[c_idx_back]] : 0.0f;

          float rho;
          if (front_valid && back_valid) {
            rho =
                0.5f * (mac.densities[c_idx_front] + mac.densities[c_idx_back]);
          } else if (front_valid) {
            rho = mac.densities[c_idx_front];
          } else {
            rho = mac.densities[c_idx_back];
          }

          float correction = (p_front - p_back) * (dt / rho);
          correction /= mac.dz;

          mac_next.w[mac.w_idx(i, j, k)] -= correction;
        } else {
          mac_next.w[mac.w_idx(i, j, k)] = 0.0f;
        }
      }
    }
  }
}

Eigen::VectorXd
Simulation3D::_build_divergence(const std::vector<int> &fluid_idx,
                                int fluid_count) {
  Eigen::VectorXd divergence(fluid_count);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      for (int k = 0; k < nz; ++k) {
        int c_idx = mac.idx(i, j, k);
        if (fluid_idx[c_idx] == -1) {
          continue;
        }

        float div = 0.0f;

        div += mac.get_cell_type(i + 1, j, k) == CellType::SOLID
                   ? 0.0f
                   : mac_next.u[mac.u_idx(i + 1, j, k)]; // right face

        div -= mac.get_cell_type(i - 1, j, k) == CellType::SOLID
                   ? 0.0f
                   : mac_next.u[mac.u_idx(i, j, k)]; // left face

        div += mac.get_cell_type(i, j + 1, k) == CellType::SOLID
                   ? 0.0f
                   : mac_next.v[mac.v_idx(i, j + 1, k)];

        div -= mac.get_cell_type(i, j - 1, k) == CellType::SOLID
                   ? 0.0f
                   : mac_next.v[mac.v_idx(i, j, k)];

        div += mac.get_cell_type(i, j, k + 1) == CellType::SOLID
                   ? 0.0f
                   : mac_next.w[mac.w_idx(i, j, k + 1)];
        div -= mac.get_cell_type(i, j, k - 1) == CellType::SOLID
                   ? 0.0f
                   : mac_next.w[mac.w_idx(i, j, k)];

        divergence[fluid_idx[c_idx]] = -div / mac.dx;

        if (std::isnan(divergence[fluid_idx[c_idx]])) {
          std::cout << "[ERROR] NaN divergence at cell (" << i << ", " << j
                    << ", " << k << ")." << std::endl;
        }

        if (std::isinf(divergence[fluid_idx[c_idx]])) {
          std::cout << "[ERROR] Inf divergence at cell (" << i << ", " << j
                    << ", " << k << ")." << std::endl;
        }
      }
    }
  }

  return divergence;
}

Eigen::SparseMatrix<double>
Simulation3D::_build_pressure_matrix(const std::vector<int> &fluid_idx,
                                     int fluid_count) {
  std::vector<Eigen::Triplet<double>> triplet_list;

  float inv_dx2 = 1.0f / (mac.dx * mac.dx);

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int c_idx = mac.idx(i, j, k);
        if (fluid_idx[c_idx] == -1) {
          continue; // not a fluid cell
        }

        int mat_row = fluid_idx[c_idx];
        int neighbor_count = 0;

        // check neighbors
        const std::vector<std::tuple<int, int, int>> neighbors = {
            {1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
            {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};

        for (const auto &[di, dj, dk] : neighbors) {
          int ni = i + di;
          int nj = j + dj;
          int nk = k + dk;

          if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz) {
            int n_idx = mac.idx(ni, nj, nk);
            if (fluid_idx[n_idx] != -1) {
              // fluid neighbor
              triplet_list.emplace_back(mat_row, fluid_idx[n_idx], -inv_dx2);
            }

            if (mac.get_cell_type(ni, nj, nk) != CellType::SOLID) {
              // non-solid neighbor
              neighbor_count++;
            }
          }
        }

        // diagonal entry
        triplet_list.emplace_back(
            mat_row, mat_row, static_cast<double>(neighbor_count) * inv_dx2);
      }
    }
  }

  Eigen::SparseMatrix<double> A(fluid_count, fluid_count);
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());
  return A;
}

void Simulation3D::_advect_cell_data(float dt) {
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) { // cell-centered
        if (mac.get_cell_type(i, j, k) == CellType::SOLID) {
          mac_next.densities[mac.idx(i, j, k)] =
              mac.densities[mac.idx(i, j, k)];
          continue;
        }

        float x = (i + 0.5f) * mac.dx;
        float y = (j + 0.5f) * mac.dy;
        float z = (k + 0.5f) * mac.dz;

        vec3d x_dest = vec3d(x, y, z);
        vec3d x_orig = IVPSolvers<vec3d, MAC3D &>::solveIVP(
            ivp_solver, x_dest, mac.current_time, dt, MAC3D::dx_vel_dt, mac);

        // if (x_orig[0] < 0.0f || x_orig[0] > nx * mac.dx || x_orig[1] < 0.0f
        // ||
        //     x_orig[1] > ny * mac.dy || x_orig[2] < 0.0f ||
        //     x_orig[2] > nz * mac.dz) {
        //   // out of bounds, set density to 0
        //   std::cout << "[WARNING] Cell data advection out of bounds at (" <<
        //   i
        //             << ", " << j << ", " << k << ")."
        //             << std::endl;
        // }

        if (mac.is_position_solid(x_orig)) {
          x_orig = mac.nonsolid_projection(x_orig, x_dest,
                                           nonsolid_projection_iterations);
        }

        mac_next.densities[mac.idx(i, j, k)] = mac.density(x_orig);
      }
    }
  }
}
