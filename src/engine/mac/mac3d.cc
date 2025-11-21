#include "mac3d.h"

MAC3d::MAC3d(int nx, int ny, int nz, float dx, float dy, float dz)
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
