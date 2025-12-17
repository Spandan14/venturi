#include "mac3d.h"
#include "utils/interpolators.h"

MAC3D::MAC3D(int nx, int ny, int nz, float dx, float dy, float dz)
    : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz) {

  u.resize((nx + 1) * ny * nz, 0.0f);
  v.resize(nx * (ny + 1) * nz, 0.0f);
  w.resize(nx * ny * (nz + 1), 0.0f);

  pressures.resize(nx * ny * nz, 0.0f);
  forces.resize(nx * ny * nz, vec3d(0.0, 0.0, 0.0));
  densities.resize(nx * ny * nz, 0.0f);

  cells.resize(nx * ny * nz);
  is_solid.resize(nx * ny * nz, false);

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
      }
    }
  }
}

float MAC3D::vel_u(vec3d pos) const {
  float x = pos[0], y = pos[1], z = pos[2];

  x = std::clamp(x, 0.0f, nx * dx);
  y = std::clamp(y, dy * 0.5f, (ny - 0.5f) * dy);
  z = std::clamp(z, dz * 0.5f, (nz - 0.5f) * dz);

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

float MAC3D::vel_v(vec3d pos) const {
  float x = pos[0], y = pos[1], z = pos[2];

  x = std::clamp(x, dx * 0.5f, (nx - 0.5f) * dx);
  y = std::clamp(y, 0.0f, ny * dy);
  z = std::clamp(z, dz * 0.5f, (nz - 0.5f) * dz);

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

float MAC3D::vel_w(vec3d pos) const {
  float x = pos[0], y = pos[1], z = pos[2];

  x = std::clamp(x, dx * 0.5f, (nx - 0.5f) * dx);
  y = std::clamp(y, dy * 0.5f, (ny - 0.5f) * dy);
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

float MAC3D::density(vec3d pos) const {
  float x = pos[0], y = pos[1], z = pos[2];

  x = std::clamp(x, 0.0f + 0.5f * dx, (nx - 0.5f) * dx);
  y = std::clamp(y, 0.0f + 0.5f * dy, (ny - 0.5f) * dy);
  z = std::clamp(z, 0.0f + 0.5f * dz, (nz - 0.5f) * dz);

  float gx = (x / dx) - 0.5f;
  float gy = (y / dy) - 0.5f;
  float gz = (z / dz) - 0.5f;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);
  int k0 = std::floor(gz);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx - 1);
  int j1 = std::min(j0 + 1, ny - 1);
  int k1 = std::min(k0 + 1, nz - 1);

  float tx = gx - i0;
  float ty = gy - j0;
  float tz = gz - k0;

  return Interpolators::trilinear(
      densities[idx(i0, j0, k0)], densities[idx(i0, j0, k1)],
      densities[idx(i0, j1, k0)], densities[idx(i0, j1, k1)],
      densities[idx(i1, j0, k0)], densities[idx(i1, j0, k1)],
      densities[idx(i1, j1, k0)], densities[idx(i1, j1, k1)], tx, ty, tz);
}

vec3d MAC3D::dx_vel_dt(vec3d x, float t, MAC3D &mac) {
  return -1.f * mac.vel(x);
}

CellType MAC3D::get_cell_type(int i, int j, int k) const {
  if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
    return CellType::SOLID; // out of bounds treated as solid
  }

  if (is_solid[idx(i, j, k)])
    return CellType::SOLID;

  return CellType::FLUID;
}

bool MAC3D::is_position_solid(vec3d pos) const {
  float x = pos[0], y = pos[1], z = pos[2];

  int i = std::clamp(static_cast<int>(x / dx), 0, nx - 1);
  int j = std::clamp(static_cast<int>(y / dy), 0, ny - 1);
  int k = std::clamp(static_cast<int>(z / dz), 0, nz - 1);

  return is_solid[idx(i, j, k)];
}

vec3d MAC3D::nonsolid_projection(vec3d in_solid, vec3d origin,
                                 int total_iter) const {
  // simple projection: move back towards origin until outside solid
  if (!is_position_solid(in_solid)) {
    return in_solid;
  }

  vec3d dir = in_solid - origin;
  float len = dir.norm();
  if (len < 1e-6f) {
    return origin; // can't determine direction, return origin
  }

  dir /= len; // normalize

  float t0 = 0.0f, t1 = 1.0f;
  vec3d temp;
  for (int iter = 0; iter < total_iter; ++iter) {
    float tm = 0.5f * (t0 + t1);
    temp = origin + dir * (len * tm);
    if (is_position_solid(temp)) {
      t1 = tm;
    } else {
      t0 = tm;
    }
  }

  return origin + dir * (len * t0);
}
