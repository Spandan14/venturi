#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

class Renderer {
public:
  Renderer() = default;
  virtual ~Renderer() = default;

  inline bool should_draw() const { return !glfwWindowShouldClose(_window); }
  virtual void render() = 0;
  inline void poll() { glfwPollEvents(); }

protected:
  GLFWwindow *_window;
};
