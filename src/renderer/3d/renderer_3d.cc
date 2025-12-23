#include "renderer_3d.h"
#include "fluxlang/utils.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include "renderer/renderer.h"

Renderer3D::Renderer3D(const MAC3D &mac, std::shared_ptr<Camera> camera,
                       int screen_width, int screen_height)
    : mac(mac), camera(camera), screen_width(screen_width),
      screen_height(screen_height), program(0), vao(0), vbo(0) {
  camera->aspect_ratio =
      static_cast<float>(screen_width) / static_cast<float>(screen_height);
  gimbal_control = std::make_shared<GimbalControl>(
      GimbalControl(this->camera,
                    vec3f(mac.nx * mac.dx / 2.0f, mac.ny * mac.dy / 2.0f,
                          mac.nz * mac.dz / 2.0f),
                    mac.nx * mac.dx * 8.0f));
  _dragging = false;
  _last_mouse_x = 0.0;
  _last_mouse_y = 0.0;

  _init_gl();
  _init_data();
  _init_handlers();
  _init_imgui();
}

Renderer3D::Renderer3D(const MAC3D &mac, int screen_width, int screen_height)
    : mac(mac), screen_width(screen_width), screen_height(screen_height),
      program(0), vao(0), vbo(0) {
  camera = std::make_shared<Camera>(Camera());
  camera->aspect_ratio =
      static_cast<float>(screen_width) / static_cast<float>(screen_height);
  gimbal_control = std::make_shared<GimbalControl>(
      GimbalControl(camera,
                    vec3f(mac.nx * mac.dx / 2.0f, mac.ny * mac.dy / 2.0f,
                          mac.nz * mac.dz / 2.0f),
                    mac.nx * mac.dx * 8.0f));
  _dragging = false;
  _last_mouse_x = 0.0;
  _last_mouse_y = 0.0;

  _init_gl();
  _init_data();
  _init_handlers();
  _init_imgui();
  _setup_imgui();
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
  glfwSetWindowUserPointer(_window, this);

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

  gimbal_control->update_camera();
  // camera->position =
  //     vec3f(mac.nx * mac.dx / 2.0f, mac.ny * mac.dy / 2.0f,
  //           -8 * std::max({mac.nx * mac.dx, mac.ny * mac.dy, mac.nz *
  //           mac.dz}));
  // camera->rotate_deg(180.f, vec3f(0.0f, 1.0f, 0.0f));
}

void Renderer3D::_init_handlers() {
  glfwSetMouseButtonCallback(_window, _mouse_button_callback);
  glfwSetCursorPosCallback(_window, _mouse_move_callback);
  glfwSetScrollCallback(_window, _mouse_scroll_callback);
  glfwSetFramebufferSizeCallback(_window, _framebuffer_resize_callback);
}

void Renderer3D::_init_imgui() {
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();

  ImGui_ImplGlfw_InitForOpenGL(_window, true);
  ImGui_ImplOpenGL3_Init("#version 410");
}

void Renderer3D::_setup_imgui() {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  ImGui::SetNextWindowSize(ImVec2(200, screen_height));
  ImGui::SetNextWindowPos(ImVec2(screen_width - 200, 0));
  ImGui::Begin("Debug");

  // ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 1, 1, 1));
  ImGui::Text("FPS: %.3f", ImGui::GetIO().Framerate);
  ImGui::Separator();

  ImGui::Text("Simulation Time: %.3f s", mac.current_time);
  ImGui::Separator();

  // ImGui::PopStyleColor();
  ImGui::End();

  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void Renderer3D::_load_uniforms() {
  glUniform3fv(glGetUniformLocation(program, "camera_pos"), 1,
               camera->position.data());
  glUniformMatrix4fv(glGetUniformLocation(program, "inv_view"), 1, GL_FALSE,
                     camera->get_view_matrix_inverse().data());
  glUniformMatrix4fv(glGetUniformLocation(program, "inv_proj"), 1, GL_FALSE,
                     camera->get_projection_matrix_inverse().data());
  glUniform1f(glGetUniformLocation(program, "aspect_ratio"),
              camera->aspect_ratio);

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

void Renderer3D::_on_mouse_button(int button, int action, int mods) {
  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if (action == GLFW_PRESS) {
      this->_dragging = true;
      glfwGetCursorPos(_window, &this->_last_mouse_x, &this->_last_mouse_y);
    } else if (action == GLFW_RELEASE) {
      this->_dragging = false;
    }
  }
}

void Renderer3D::_mouse_button_callback(GLFWwindow *window, int button,
                                        int action, int mods) {
  auto *renderer = static_cast<Renderer3D *>(glfwGetWindowUserPointer(window));

  if (!renderer) {
    return;
  }

  renderer->_on_mouse_button(button, action, mods);
}

void Renderer3D::_on_mouse_move(double x_pos, double y_pos) {
  if (!_dragging) {
    return;
  }

  double delta_x = x_pos - this->_last_mouse_x;
  double delta_y = y_pos - this->_last_mouse_y;

  this->gimbal_control->handle_mouse_move(static_cast<float>(delta_x),
                                          static_cast<float>(delta_y));

  this->_last_mouse_x = x_pos;
  this->_last_mouse_y = y_pos;
}

void Renderer3D::_mouse_move_callback(GLFWwindow *window, double x_pos,
                                      double y_pos) {
  auto *renderer = static_cast<Renderer3D *>(glfwGetWindowUserPointer(window));

  if (!renderer) {
    return;
  }

  renderer->_on_mouse_move(x_pos, y_pos);
}

void Renderer3D::_on_mouse_scroll(double x_offset, double y_offset) {
  this->gimbal_control->handle_mouse_scroll(static_cast<float>(y_offset));
}

void Renderer3D::_mouse_scroll_callback(GLFWwindow *window, double x_offset,
                                        double y_offset) {
  auto *renderer = static_cast<Renderer3D *>(glfwGetWindowUserPointer(window));

  if (!renderer) {
    return;
  }

  renderer->_on_mouse_scroll(x_offset, y_offset);
}

void Renderer3D::_on_framebuffer_resize(int width, int height) {
  this->screen_width = width;
  this->screen_height = height;
  glViewport(0, 0, width, height);

  camera->aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
}

void Renderer3D::_framebuffer_resize_callback(GLFWwindow *window, int width,
                                              int height) {
  auto *renderer = static_cast<Renderer3D *>(glfwGetWindowUserPointer(window));

  if (!renderer) {
    return;
  }

  renderer->_on_framebuffer_resize(width, height);
}

Renderer3D::~Renderer3D() {
  // imgui shutdown
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

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
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glDisable(GL_BLEND);
  _load_uniforms();
  _load_sim_data();

  _draw_sim();

  glUseProgram(0);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  _setup_imgui();

  glfwSwapBuffers(_window);
}

void Renderer3D::_draw_sim() {
  glBindVertexArray(vao);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}
