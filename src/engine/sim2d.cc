#include "sim2d.h"
#include <utils/physical_consts.h>

Simulation2D::Simulation2D(int nx, int ny, float dx, float dy)
    : mac(nx, ny, dx, dy), mac_next(nx, ny, dx, dy), nx(nx), ny(ny), dx(dx),
      dy(dy) {}

void Simulation2D::step(float dt) {
  _advect_velocities(dt);
  _advect_cell_data(dt);

  mac_next.current_time += dt;

  _apply_forces(dt);

  // everything has now ended up in mac_next, so we swap
  std::swap(mac, mac_next);
}

void Simulation2D::initialize_vel_u(
    const std::function<float(int, int)> &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) {
      mac.u[mac.u_idx(i, j)] = initializer(i, j);
    }
  }
}

void Simulation2D::initialize_vel_v(
    const std::function<float(int, int)> &initializer) {
  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.v[mac.v_idx(i, j)] = initializer(i, j);
    }
  }
}

void Simulation2D::initialize_density(
    const std::function<float(int, int)> &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.densities[mac.idx(i, j)] = initializer(i, j);
    }
  }
}

void Simulation2D::initialize_forces(
    const std::function<vec2d(int, int)> &initializer) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      mac.forces[mac.idx(i, j)] = initializer(i, j);
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
          solver, x_dest, mac.current_time, dt, MAC2D::dx_vel_dt, mac);

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
          solver, x_dest, mac.current_time, dt, MAC2D::dx_vel_dt, mac);

      // interpolate the velocity at the original position
      mac_next.v[mac.v_idx(i, j)] = mac.vel_v(x_orig);
    }
  }
}

void Simulation2D::_advect_cell_data(float dt) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) { // cell-centered
      float x = (i + 0.5f) * mac.dx;
      float y = (j + 0.5f) * mac.dy;

      vec2d x_dest = vec2d(x, y);

      // solve for original position
      vec2d x_orig = IVPSolvers<vec2d, MAC2D &>::solveIVP(
          solver, x_dest, mac.current_time, dt, MAC2D::dx_vel_dt, mac);

      // interpolate the cell data at the original position
      // for now, just density
      mac_next.densities[mac.idx(i, j)] = mac.density(x_orig);
    }
  }
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

      mac_next.u[mac.u_idx(i, j)] += force * dt;
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

      mac_next.v[mac.v_idx(i, j)] += force * dt;
    }
  }
}
