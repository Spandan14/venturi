#pragma once

#include <engine/mac/mac2d.h>

class Renderer2D {
public:
  Renderer2D(MAC2D &mac, int screen_width, int screen_height);
  ~Renderer2D();

  void render();

private:
  MAC2D &mac;

  int screen_width, screen_height;
  GLuint program;
  GLuint vao, vbo;

  GLuint compile_slang_shader
}
