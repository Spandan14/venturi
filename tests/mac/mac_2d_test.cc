#include "mac/mac2d.h"
#include <iostream>

int test_one() {
  int nx = 4, ny = 4, dx = 1, dy = 1;
  MAC2D mac(nx, ny, dx, dy);

  // affine velocity field: u = x + 2y, v = -x + y
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i <= nx; ++i) { // u-grid
      float x = i * mac.dx;
      float y = (j + 0.5f) * mac.dy;
      mac.u[mac.u_idx(i, j)] = x + 2 * y;
    }
  }
  for (int j = 0; j <= ny; ++j) {
    for (int i = 0; i < nx; ++i) { // v-grid
      float x = (i + 0.5f) * mac.dx;
      float y = j * mac.dy;
      mac.v[mac.v_idx(i, j)] = -x + y;
    }
  }

  // test points
  for (float y = 0.1f; y < ny; y += 0.03f) {
    for (float x = 0.1f; x < nx; x += 0.03f) {

      float u_exact = x + 2 * y;
      float v_exact = -x + y;

      float u_interp = mac.u_vel(x, y);
      float v_interp = mac.v_vel(x, y);

      assert(std::abs(u_interp - u_exact) < 1e-6f);
      assert(std::abs(v_interp - v_exact) < 1e-6f);
    }
  }

  std::cout << "All affine interpolation tests passed!" << std::endl;
  return 0;
}

int main() {
  test_one();
  return 0;
}
