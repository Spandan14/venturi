#pragma once

#include "renderer/renderer.h"
#include <Eigen/Core>
#include <engine/mac/mac2d.h>

typedef Eigen::Vector3d vec3d;

class Renderer2D : public Renderer {
public:
  Renderer2D(const MAC2D &mac, int screen_width, int screen_height);
  ~Renderer2D();

  void render() override;

private:
  const MAC2D &mac;

  int screen_width, screen_height;
  GLuint program;
  GLuint vao, vbo;

  void _draw_sim();

  void _draw_outline(float x0, float y0, float x1, float y1, float r, float g,
                     float b);
  void _draw_quad(float x0, float y0, float x1, float y1, float r, float g,
                  float b);
};
