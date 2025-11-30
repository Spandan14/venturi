#pragma once

class Simulation {
public:
  Simulation() = default;
  virtual ~Simulation() = default;

  virtual void step(float dt) = 0;
};
