#include "renderer_3d.h"
#include "fluxlang/utils.h"
#include "renderer/renderer.h"
#include "utils/interpolators.h"
#include <iostream>

Renderer3D::Renderer3D(const MAC3D &mac, const Camera &camera, int screen_width,
                       int screen_height)
    : mac(mac), camera(camera), screen_width(screen_width),
      screen_height(screen_height), program(0), vao(0), vbo(0) {
  _init_gl();
  _init_data();
}

Renderer3D::Renderer3D(const MAC3D &mac, int screen_width, int screen_height)
    : mac(mac), screen_width(screen_width), screen_height(screen_height),
      program(0), vao(0), vbo(0) {
  camera = Camera();
  _init_gl();
  _init_data();
}

void Renderer3D::_init_gl() {
  if (!glfwInit()) {
    throw std::runtime_error("Failed to initialize GLFW");
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  _window = glfwCreateWindow(screen_width, screen_height, "Venturi 1", nullptr,
                             nullptr);
  if (!_window) {
    glfwTerminate();
    throw std::runtime_error("Failed to create GLFW window");
  }

  glfwMakeContextCurrent(_window);

  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    throw std::runtime_error("Failed to initialize GLEW");
  }

  auto vs_source = load_text_file(std::string(SHADERS_3D) + "3d.vert");
  auto fs_source = load_text_file(std::string(SHADERS_3D) + "3d.frag");

  GLuint vertex_shader = compile_shader(vs_source, GL_VERTEX_SHADER);
  GLuint fragment_shader = compile_shader(fs_source, GL_FRAGMENT_SHADER);

  program = glCreateProgram();
  glAttachShader(program, vertex_shader);
  glAttachShader(program, fragment_shader);
  glLinkProgram(program);

  GLint success;
  glGetProgramiv(program, GL_LINK_STATUS, &success);
  if (!success) {
    char info_log[512];
    glGetProgramInfoLog(program, 512, nullptr, info_log);
    std::string error_msg =
        std::string("Shader program linking failed: ") + info_log;
    throw std::runtime_error(error_msg);
  }

  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);

  if (glewInit() != GLEW_OK) {
    throw std::runtime_error("Failed to initialize GLEW");
  }

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
}

void Renderer3D::_init_data() {
  glGenVertexArrays(1, &vao);
  glGenBuffers(1, &vbo);

  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  float quad[] = {-1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, 1.f};
  glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);

  glGenTextures(1, &density_texture);

  camera.position =
      vec3f(mac.nx * mac.dx / 2.0f, mac.ny * mac.dy / 2.0f,
            -5 * std::max({mac.nx * mac.dx, mac.ny * mac.dy, mac.nz * mac.dz}));
  camera.rotate_deg(180.f, vec3f(0.0f, 1.0f, 0.0f));
}

void Renderer3D::_load_uniforms() {
  glUniform3fv(glGetUniformLocation(program, "camera_pos"), 1,
               camera.position.data());
  glUniformMatrix4fv(glGetUniformLocation(program, "inv_view"), 1, GL_FALSE,
                     camera.get_view_matrix_inverse().data());
  glUniformMatrix4fv(glGetUniformLocation(program, "inv_proj"), 1, GL_FALSE,
                     camera.get_projection_matrix_inverse().data());

  glUniform3f(glGetUniformLocation(program, "volume_min"), 0.0f, 0.0f, 0.0f);
  glUniform3f(glGetUniformLocation(program, "volume_max"), mac.nx * mac.dx,
              mac.ny * mac.dy, mac.nz * mac.dz);
  glUniform3i(glGetUniformLocation(program, "grid_size"), mac.nx, mac.ny,
              mac.nz);

  glUniform1f(glGetUniformLocation(program, "step_size"), _rm_step_size);
  glUniform1f(glGetUniformLocation(program, "absorption"), _rm_absorption);
  glUniform1i(glGetUniformLocation(program, "max_steps"), _rm_max_steps);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_3D, density_texture);
  glUniform1i(glGetUniformLocation(program, "density_tex"), 0);
}

void Renderer3D::_load_sim_data() {
  glBindTexture(GL_TEXTURE_3D, density_texture);

  glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, mac.nx, mac.ny, mac.nz, 0, GL_RED,
               GL_FLOAT, mac.densities.data());

  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
}

Renderer3D::~Renderer3D() {
  glDeleteProgram(program);
  glfwDestroyWindow(_window);
  glfwTerminate();
}

GLuint Renderer3D::compile_shader(const std::string &source,
                                  GLenum shader_type) {
  GLuint shader = glCreateShader(shader_type);
  const char *source_cstr = source.c_str();
  glShaderSource(shader, 1, &source_cstr, nullptr);
  glCompileShader(shader);

  GLint success;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  if (!success) {
    char info_log[512];
    glGetShaderInfoLog(shader, 512, nullptr, info_log);
    std::string error_msg =
        std::string("Shader compilation failed: ") + info_log;
    throw std::runtime_error(error_msg);
  }

  return shader;
}

void Renderer3D::render() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glUseProgram(program);
  _load_uniforms();
  _load_sim_data();

  // move camera in a circle around the center of the volume, in a x-z axis ring
  float time = glfwGetTime();
  float radius =
      8 * std::max({mac.nx * mac.dx, mac.ny * mac.dy, mac.nz * mac.dz});
  camera.position = vec3f(mac.nx * mac.dx / 2.0f, mac.ny * mac.dy / 2.0f,
                          mac.nz * mac.dz / 2.0f) +
                    vec3f(radius * std::cos(time * 0.2f), 0.0f,
                          radius * std::sin(time * 0.2f));
  camera.look_at(vec3f(mac.nx * mac.dx / 2.0f, mac.ny * mac.dy / 2.0f,
                       mac.nz * mac.dz / 2.0f));

  _draw_sim();

  glfwSwapBuffers(_window);
}

void Renderer3D::_draw_sim() {
  glBindVertexArray(vao);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}
