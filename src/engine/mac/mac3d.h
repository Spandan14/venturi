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

class MAC3d {

public:
  MAC3d(int nx, int ny, int nz, float dx, float dy, float dz);
  ~MAC3d();

  int nx, ny, nz;
  float dx, dy, dz;

private:
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
