#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <engine/mac/mac2d.h>

class Renderer2D {
public:
  Renderer2D(const MAC2D &mac, int screen_width, int screen_height);
  ~Renderer2D();

  inline bool should_draw() const { return !glfwWindowShouldClose(_window); }
  void render();
  inline void poll() { glfwPollEvents(); }

private:
  const MAC2D &mac;

  int screen_width, screen_height;
  GLuint program;
  GLuint vao, vbo;

  GLFWwindow *_window;

  void _draw_sim();

  void _draw_outline(float x0, float y0, float x1, float y1);
  void _draw_quad(float x0, float y0, float x1, float y1, float r, float g,
                  float b);
};
