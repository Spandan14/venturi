#include "mac2d.h"

MAC2d::MAC2d(int nx, int ny, float dx, float dy)
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
