#include "sim2d.h"

Simulation2D::Simulation2D(int nx, int ny, float dx, float dy)
    : mac(nx, ny, dx, dy), nx(nx), ny(ny), dx(dx), dy(dy) {}
