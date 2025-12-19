#pragma once

#include "renderer/3d/camera.h"
#include "renderer/renderer.h"
#include <Eigen/Core>
#include <engine/mac/mac3d.h>

typedef Eigen::Vector3d vec3d;

class Renderer3D : public Renderer {
public:
  Renderer3D(const MAC3D &mac, const Camera &camera, int screen_width,
             int screen_height);
  Renderer3D(const MAC3D &mac, int screen_width, int screen_height);
  ~Renderer3D();

  void render() override;

private:
  const MAC3D &mac;
  Camera camera;

  int screen_width, screen_height;
  GLuint program;
  GLuint vao, vbo;
  GLuint density_texture;

  float _rm_step_size = 0.1f;
  float _rm_absorption = 1.0f;
  int _rm_max_steps = 256;

  void _init_gl();
  void _init_data();

  void _load_uniforms();
  void _load_sim_data();

  GLuint compile_shader(const std::string &source, GLenum shader_type);

  void _draw_sim();
};
