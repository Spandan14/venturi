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

struct CellData3D {
  float pressure;
  vec3d force;
  float density;
};

class MAC3D {

public:
  MAC3D(int nx, int ny, int nz, float dx, float dy, float dz);
  ~MAC3D() = default;

  int nx, ny, nz;
  float dx, dy, dz;

  float u_vel(float x, float y, float z) const;
  float v_vel(float x, float y, float z) const;
  float w_vel(float x, float y, float z) const;

  inline int u_idx(int i, int j, int k) const {
    return i + j * (nx + 1) + k * (nx + 1) * ny;
  }
  inline int v_idx(int i, int j, int k) const {
    return i + j * nx + k * nx * (ny + 1);
  }
  inline int w_idx(int i, int j, int k) const {
    return i + j * nx + k * nx * ny;
  }

  std::vector<Cell3D> cells;

  std::vector<float> u; // x-face velocities
                        // size: (nx + 1) * ny * nz

  std::vector<float> v; // y-face velocities
                        // size: nx * (ny + 1) * nz

  std::vector<float> w; // z-face velocities
                        // size: nx * ny * (nz + 1)

  std::vector<CellData3D> cell_data; // cell-centered data
                                     // size: nx * ny * nz
};
