#include "mac2d.h"
#include "utils/interpolators.h"

MAC2D::MAC2D(int nx, int ny, float dx, float dy)
    : nx(nx), ny(ny), dx(dx), dy(dy) {

  u.resize((nx + 1) * ny, 0.0f);
  v.resize(nx * (ny + 1), 0.0f);

  pressures.resize(nx * ny, 0.0f);
  forces.resize(nx * ny, vec2d(0.0, 0.0));
  densities.resize(nx * ny, 0.0f);

  cells.resize(nx * ny);
  is_solid.resize(nx * ny, false);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int c_idx = i + j * nx;
      Cell2D &cell = cells[c_idx];
      cell.c_idx = c_idx;

      cell.u_lo_idx = i + j * (nx + 1);
      cell.u_hi_idx = (i + 1) + j * (nx + 1);

      cell.v_lo_idx = i + j * nx;
      cell.v_hi_idx = i + (j + 1) * nx;
    }
  }
}

float MAC2D::vel_u(vec2d pos) const {
  float x = pos[0], y = pos[1];

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

float MAC2D::vel_v(vec2d pos) const {
  float x = pos[0], y = pos[1];

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

float MAC2D::density(vec2d pos) const {
  float x = pos[0], y = pos[1];

  x = std::clamp(x, 0.0f + 0.5f * dx, (nx - 0.5f) * dx);
  y = std::clamp(y, 0.0f + 0.5f * dy, (ny - 0.5f) * dy);

  float gx = (x / dx) - 0.5f;
  float gy = (y / dy) - 0.5f;

  int i0 = std::floor(gx);
  int j0 = std::floor(gy);

  // WARN: does this work?
  int i1 = std::min(i0 + 1, nx - 1);
  int j1 = std::min(j0 + 1, ny - 1);

  float tx = gx - i0;
  float ty = gy - j0;

  return Interpolators::bilinear(densities[idx(i0, j0)], densities[idx(i0, j1)],
                                 densities[idx(i1, j0)], densities[idx(i1, j1)],
                                 tx, ty);
}

vec2d MAC2D::dx_vel_dt(vec2d x, float t, MAC2D &mac) {
  return -1.f * mac.vel(x);
}

CellType MAC2D::get_cell_type(int i, int j) const {
  if (i < 0 || i >= nx || j < 0 || j >= ny)
    return CellType::SOLID; // out of bounds treated as solid

  // for now, everything is fluid
  if (is_solid[idx(i, j)])
    return CellType::SOLID;

  // if (densities[idx(i, j)] > 0.0f)
  return CellType::FLUID;

  // return CellType::AIR;
}

bool MAC2D::is_position_solid(vec2d pos) const {
  float x = pos[0], y = pos[1];

  int i = std::clamp(int(x / dx), 0, nx - 1);
  int j = std::clamp(int(y / dy), 0, ny - 1);

  return is_solid[idx(i, j)];
}

vec2d MAC2D::nonsolid_projection(vec2d in_solid, vec2d origin,
                                 int total_iter) const {
  // simple projection: move back towards origin until outside solid
  if (!is_position_solid(in_solid)) {
    return in_solid;
  }

  vec2d dir = in_solid - origin;
  float len = dir.norm();
  if (len < 1e-6f) {
    return origin; // can't determine direction, return origin
  }

  dir /= len; // normalize

  float t0 = 0.0f, t1 = 1.0f;
  vec2d temp;
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
