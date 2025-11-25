#pragma once

#include "peglib.h"
#include <engine/sim2d.h>

class Parser2D {
public:
  Parser2D() = default;
  ~Parser2D() = default;

  static void parse(const char *filename, Simulation2D &sim);
};
