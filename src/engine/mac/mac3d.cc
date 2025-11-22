#include "mac3d.h"
#include "utils/interpolators.h"

MAC3D::MAC3D(int nx, int ny, int nz, float dx, float dy, float dz)
    : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz) {

  cells.resize(nx * ny * nz);
  cell_data.resize(nx * ny * nz);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int c_idx = i + j * nx + k * nx * ny;
        Cell3D &cell = cells[c_idx];
        cell.c_idx = c_idx;

        cell.u_lo_idx = i + j * (nx + 1) + k * (nx + 1) * ny;
        cell.u_hi_idx = (i + 1) + j * (nx + 1) + k * (nx + 1) * ny;

        cell.v_lo_idx = i + j * nx + k * nx * (ny + 1);
        cell.v_hi_idx = i + (j + 1) * nx + k * nx * (ny + 1);

        cell.w_lo_idx = i + j * nx + k * nx * ny;
        cell.w_hi_idx = i + j * nx + (k + 1) * nx * ny;

        CellData3D &data = cell_data[c_idx];
        data.pressure = 0.0f;
        data.force = vec3d(0.0, 0.0, 0.0);
        data.density = 0.0f;
      }
    }
  }

  u.resize((nx + 1) * ny * nz, 0.0f);
  v.resize(nx * (ny + 1) * nz, 0.0f);
  w.resize(nx * ny * (nz + 1), 0.0f);
}

float MAC3D::u_vel(float x, float y, float z) const {
  x = std::clamp(x, 0.0f, nx * dx);
  y = std::clamp(y, dy * 0.5f,
                 (ny - 0.5f) * dy); // because u velocities are at x-faces
  z = std::clamp(z, dz * 0.5f,
                 (nz - 0.5f) * dz); // because u velocities are at x-faces

  float gx = x / dx;
  float gy = (y / dy) - 0.5f;
  float gz = (z / dz) - 0.5f;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);
  int k0 = std::floor(gz);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx);
  int j1 = std::min(j0 + 1, ny - 1);
  int k1 = std::min(k0 + 1, nz - 1);

  float tx = gx - i0;
  float ty = gy - j0;
  float tz = gz - k0;

  return Interpolators::trilinear(
      u[u_idx(i0, j0, k0)], u[u_idx(i0, j0, k1)], u[u_idx(i0, j1, k0)],
      u[u_idx(i0, j1, k1)], u[u_idx(i1, j0, k0)], u[u_idx(i1, j0, k1)],
      u[u_idx(i1, j1, k0)], u[u_idx(i1, j1, k1)], tx, ty, tz);
}

float MAC3D::v_vel(float x, float y, float z) const {
  x = std::clamp(x, dx * 0.5f,
                 (nx - 0.5f) * dx); // because v velocities are at y-faces
  y = std::clamp(y, 0.0f, ny * dy);
  z = std::clamp(z, dz * 0.5f,
                 (nz - 0.5f) * dz); // because v velocities are at y-faces

  float gx = (x / dx) - 0.5f;
  float gy = y / dy;
  float gz = (z / dz) - 0.5f;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);
  int k0 = std::floor(gz);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx - 1);
  int j1 = std::min(j0 + 1, ny);
  int k1 = std::min(k0 + 1, nz - 1);

  float tx = gx - i0;
  float ty = gy - j0;
  float tz = gz - k0;

  return Interpolators::trilinear(
      v[v_idx(i0, j0, k0)], v[v_idx(i0, j0, k1)], v[v_idx(i0, j1, k0)],
      v[v_idx(i0, j1, k1)], v[v_idx(i1, j0, k0)], v[v_idx(i1, j0, k1)],
      v[v_idx(i1, j1, k0)], v[v_idx(i1, j1, k1)], tx, ty, tz);
}

float MAC3D::w_vel(float x, float y, float z) const {
  x = std::clamp(x, dx * 0.5f,
                 (nx - 0.5f) * dx); // because w velocities are at z-faces
  y = std::clamp(y, dy * 0.5f,
                 (ny - 0.5f) * dy); // because w velocities are at z-faces
  z = std::clamp(z, 0.0f, nz * dz);

  float gx = (x / dx) - 0.5f;
  float gy = (y / dy) - 0.5f;
  float gz = z / dz;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);
  int k0 = std::floor(gz);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx - 1);
  int j1 = std::min(j0 + 1, ny - 1);
  int k1 = std::min(k0 + 1, nz);

  float tx = gx - i0;
  float ty = gy - j0;
  float tz = gz - k0;

  return Interpolators::trilinear(
      w[w_idx(i0, j0, k0)], w[w_idx(i0, j0, k1)], w[w_idx(i0, j1, k0)],
      w[w_idx(i0, j1, k1)], w[w_idx(i1, j0, k0)], w[w_idx(i1, j0, k1)],
      w[w_idx(i1, j1, k0)], w[w_idx(i1, j1, k1)], tx, ty, tz);
}
