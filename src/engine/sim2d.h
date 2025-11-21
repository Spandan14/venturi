#pragma once
#include "mac/mac2d.h"

class Simulation2D {
public:
  Simulation2D(int nx, int ny, float dx, float dy);
  ~Simulation2D();

  void step(float dt);

private:
  MAC2D mac;
  int nx, ny;
  float dx, dy;
};
