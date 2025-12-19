#pragma once
#include <Eigen/Core>
#include <engine/mac/mac.h>
#include <vector>

typedef Eigen::Vector2d vec2d;

struct Cell2D {
  int c_idx;

  int u_lo_idx;
  int u_hi_idx;

  int v_lo_idx;
  int v_hi_idx;
};

class MAC2D {
public:
  MAC2D(int nx, int ny, float dx, float dy);
  ~MAC2D() = default;

  int nx, ny;
  float dx, dy;

  float vel_u(vec2d pos) const;
  float vel_v(vec2d pos) const;

  inline vec2d vel(vec2d pos) const { return vec2d(vel_u(pos), vel_v(pos)); }

  float density(vec2d pos) const;

  inline int u_idx(int i, int j) const { return i + j * (nx + 1); }
  inline int v_idx(int i, int j) const { return i + j * nx; }
  inline int idx(int i, int j) const { return i + j * nx; }

  float current_time = 0.0f;
  CellType get_cell_type(int i, int j) const;
  bool is_position_solid(vec2d pos) const;

  vec2d nonsolid_projection(vec2d in_solid, vec2d origin, int total_iter) const;

  std::vector<Cell2D> cells;
  std::vector<bool> is_solid; // cell-centered solid flags
                              // size: nx * ny

  std::vector<float> u; // x-face velocities
                        // size: (nx + 1) * ny

  std::vector<float> v; // y-face velocities
                        // size: nx * (ny + 1)

  std::vector<float> pressures; // cell-centered pressures
                                // size: nx * ny

  std::vector<vec2d> forces; // cell-centered forces
                             // size: nx * ny

  std::vector<float> densities; // cell-centered densities
                                // size: nx * ny

  static vec2d dx_vel_dt(vec2d x, float t, MAC2D &mac);
};
