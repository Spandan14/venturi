#pragma once
#include <Eigen/Core>
#include <vector>

typedef Eigen::Vector2d vec2d;

struct Cell2D {
  int c_idx;

  int u_lo_idx;
  int u_hi_idx;

  int v_lo_idx;
  int v_hi_idx;
};

struct CellData2D {
  float pressure;
  vec2d force;
  float density;
};

class MAC2d {
public:
  MAC2d(int nx, int ny, float dx, float dy);
  ~MAC2d();

  int nx, ny;
  float dx, dy;

private:
  std::vector<Cell2D> cells;

  std::vector<float> u; // x-face velocities
                        // size: (nx + 1) * ny

  std::vector<float> v; // y-face velocities
                        // size: nx * (ny + 1)

  std::vector<CellData2D> cell_data; // cell-centered data
                                     // size: nx * ny
};
