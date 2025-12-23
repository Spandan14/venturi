#pragma once

#include "renderer/3d/camera/camera.h"
#include "renderer/3d/camera/gimbal_control.h"
#include "renderer/renderer.h"
#include <Eigen/Core>
#include <engine/mac/mac3d.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <memory>

typedef Eigen::Vector3d vec3d;

class Renderer3D : public Renderer {
public:
  Renderer3D(const MAC3D &mac, std::shared_ptr<Camera> camera, int screen_width,
             int screen_height);
  Renderer3D(const MAC3D &mac, int screen_width, int screen_height);
  ~Renderer3D();

  void render() override;

private:
  const MAC3D &mac;
  std::shared_ptr<Camera> camera;
  std::shared_ptr<GimbalControl> gimbal_control;

  int screen_width, screen_height;
  GLuint program;
  GLuint vao, vbo;
  GLuint density_texture;

  float _rm_step_size = 0.3f;
  float _rm_absorption = 1.0f;
  int _rm_max_steps = 256;

  void _init_gl();
  void _init_data();
  void _init_handlers();
  void _init_imgui();

  void _load_uniforms();
  void _load_sim_data();

  void _setup_imgui();

  bool _dragging;
  double _last_mouse_x, _last_mouse_y;
  void _on_mouse_button(int button, int action, int mods);
  void _on_mouse_move(double x_pos, double y_pos);
  void _on_mouse_scroll(double x_offset, double y_offset);
  void _on_framebuffer_resize(int width, int height);

  static void _mouse_button_callback(GLFWwindow *window, int button, int action,
                                     int mods);
  static void _mouse_move_callback(GLFWwindow *window, double x_pos,
                                   double y_pos);
  static void _mouse_scroll_callback(GLFWwindow *window, double x_offset,
                                     double y_offset);
  static void _framebuffer_resize_callback(GLFWwindow *window, int width,
                                           int height);

  GLuint compile_shader(const std::string &source, GLenum shader_type);

  void _draw_sim();

  void _handle_input();
};
