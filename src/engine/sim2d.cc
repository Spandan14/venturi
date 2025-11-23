#include "sim2d.h"

Simulation2D::Simulation2D(int nx, int ny, float dx, float dy)
    : mac(nx, ny, dx, dy), mac_next(nx, ny, dx, dy), nx(nx), ny(ny), dx(dx),
      dy(dy) {}

void Simulation2D::step(float dt) {
  _advect_velocities(dt);

  // everything has now ended up in mac_next, so we swap
  std::swap(mac, mac_next);
}

// advects velocities using the semi-Lagrangian methods
// looks at every destination velocity grid point by axis, solves for the
// original, and advects new velocity
void Simulation2D::_advect_velocities(float dt) { _advect_u(dt); }

void Simulation2D::_advect_u(float dt) {
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) { // u-grid
      float x = i * mac.dx;
      float y = (j + 0.5f) * mac.dy;

      // trace back to find the original position
      float x_orig = x - mac.u_vel(x, y) * dt;
      float y_orig = y - mac.v_vel(x, y) * dt;

      // clamp to domain
      x_orig = std::clamp(x_orig, 0.0f, nx * dx);
      y_orig = std::clamp(y_orig, 0.0f, ny * dy);

      // interpolate the velocity at the original position
      mac_next.u[mac.u_idx(i, j)] = mac.u_vel(x_orig, y_orig);
    }
  }
}
