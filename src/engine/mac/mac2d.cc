#include "mac2d.h"
#include "utils/interpolators.h"

MAC2D::MAC2D(int nx, int ny, float dx, float dy)
    : nx(nx), ny(ny), dx(dx), dy(dy) {

  cells.resize(nx * ny);
  cell_data.resize(nx * ny);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int c_idx = i + j * nx;
      Cell2D &cell = cells[c_idx];
      cell.c_idx = c_idx;

      cell.u_lo_idx = i + j * (nx + 1);
      cell.u_hi_idx = (i + 1) + j * (nx + 1);

      cell.v_lo_idx = i + j * nx;
      cell.v_hi_idx = i + (j + 1) * nx;

      CellData2D &data = cell_data[c_idx];
      data.pressure = 0.0f;
      data.force = vec2d(0.0, 0.0);
      data.density = 0.0f;
    }
  }

  u.resize((nx + 1) * ny, 0.0f);
  v.resize(nx * (ny + 1), 0.0f);
}

float MAC2D::u_vel(float x, float y) const {
  x = std::clamp(x, 0.0f, nx * dx);
  y = std::clamp(y, dy * 0.5f,
                 (ny - 0.5f) * dy); // because u velocities are at x-faces

  float gx = x / dx;
  float gy = (y / dy) - 0.5f;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx);
  int j1 = std::min(j0 + 1, ny - 1);

  float tx = gx - i0;
  float ty = gy - j0;

  return Interpolators::bilinear(u[u_idx(i0, j0)], u[u_idx(i0, j1)],
                                 u[u_idx(i1, j0)], u[u_idx(i1, j1)], tx, ty);
}

float MAC2D::v_vel(float x, float y) const {
  x = std::clamp(x, dx * 0.5f,
                 (nx - 0.5f) * dx); // because v velocities are at y-faces
  y = std::clamp(y, 0.0f, ny * dy);

  float gx = (x / dx) - 0.5f;
  float gy = y / dy;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx - 1);
  int j1 = std::min(j0 + 1, ny);

  float tx = gx - i0;
  float ty = gy - j0;

  return Interpolators::bilinear(v[v_idx(i0, j0)], v[v_idx(i0, j1)],
                                 v[v_idx(i1, j0)], v[v_idx(i1, j1)], tx, ty);
}
