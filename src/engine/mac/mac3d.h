#pragma once
#include <Eigen/Core>
#include <vector>

typedef Eigen::Vector3d vec3d;

struct Cell3D {
  int c_idx;

  int u_lo_idx;
  int u_hi_idx;

  int v_lo_idx;
  int v_hi_idx;

  int w_lo_idx;
  int w_hi_idx;
};

enum class CellType { FLUID, SOLID, AIR };

class MAC3D {
public:
  MAC3D(int nx, int ny, int nz, float dx, float dy, float dz);
  ~MAC3D() = default;

  int nx, ny, nz;
  float dx, dy, dz;

  float vel_u(vec3d pos) const;
  float vel_v(vec3d pos) const;
  float vel_w(vec3d pos) const;

  inline vec3d vel(vec3d pos) const {
    return vec3d(vel_u(pos), vel_v(pos), vel_w(pos));
  }

  float density(vec3d pos) const;

  inline int u_idx(int i, int j, int k) const {
    return i + j * (nx + 1) + k * (nx + 1) * ny;
  }
  inline int v_idx(int i, int j, int k) const {
    return i + j * nx + k * nx * (ny + 1);
  }
  inline int w_idx(int i, int j, int k) const {
    return i + j * nx + k * nx * ny;
  }
  inline int idx(int i, int j, int k) const { return i + j * nx + k * nx * ny; }

  float current_time = 0.0f;
  CellType get_cell_type(int i, int j, int k) const;
  bool is_position_solid(vec3d pos) const;

  vec3d nonsolid_projection(vec3d in_solid, vec3d origin, int total_iter) const;

  std::vector<Cell3D> cells;
  std::vector<bool> is_solid; // cell-centered solid flags
                              // size: nx * ny * nz

  std::vector<float> u; // x-face velocities
                        // size: (nx + 1) * ny * nz

  std::vector<float> v; // y-face velocities
                        // size: nx * (ny + 1) * nz

  std::vector<float> w; // z-face velocities
                        // size: nx * ny * (nz + 1)

  std::vector<float> pressures; // cell-centered pressures
                                // size: nx * ny * nz

  std::vector<vec3d> forces; // cell-centered forces
                             // size: nx * ny * nz

  std::vector<float> densities; // cell-centered densities
                                // size: nx * ny * nz

  static vec3d dx_vel_dt(vec3d x, float t, MAC3D &mac);
};
