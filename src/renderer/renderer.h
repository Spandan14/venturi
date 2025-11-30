#pragma once

class Renderer {
public:
  Renderer() = default;
  virtual ~Renderer() = default;

  virtual void render() = 0;
};
