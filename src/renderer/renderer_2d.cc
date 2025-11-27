#include "renderer_2d.h"
#include "utils/interpolators.h"
#include <iostream>

Renderer2D::Renderer2D(const MAC2D &mac, int screen_width, int screen_height)
    : mac(mac), screen_width(screen_width), screen_height(screen_height),
      program(0), vao(0), vbo(0) {

  if (!glfwInit()) {
    throw std::runtime_error("Failed to initialize GLFW");
  }

  _window = glfwCreateWindow(screen_width, screen_height, "Venturi 1", nullptr,
                             nullptr);
  if (!_window) {
    glfwTerminate();
    throw std::runtime_error("Failed to create GLFW window");
  }

  glfwMakeContextCurrent(_window);

  if (glewInit() != GLEW_OK) {
    throw std::runtime_error("Failed to initialize GLEW");
  }

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
}

Renderer2D::~Renderer2D() {
  glfwDestroyWindow(_window);
  glfwTerminate();
}

void Renderer2D::render() {
  glClear(GL_COLOR_BUFFER_BIT);

  _draw_sim();

  glfwSwapBuffers(_window);
}

void Renderer2D::_draw_sim() {
  float W = mac.nx * mac.dx;
  float H = mac.ny * mac.dy;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, W, 0, H, -1, 1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  for (int j = 0; j < mac.ny; j++) {
    for (int i = 0; i < mac.nx; i++) {

      float d = mac.densities[mac.idx(i, j)];
      d = std::clamp(d, 0.0f, 1.0f);

      float x0 = i * mac.dx;
      float y0 = j * mac.dy;
      float x1 = x0 + mac.dx;
      float y1 = y0 + mac.dy;

      Cell2D cell = mac.cells[mac.idx(i, j)];
      vec2d cell_vel =
          vec2d(0.5f * (mac.u[cell.u_lo_idx] + mac.u[cell.u_hi_idx]),
                0.5f * (mac.v[cell.v_lo_idx] + mac.v[cell.v_hi_idx]));

      // std::cout << "Cell (" << i << ", " << j << ") velocity: (" <<
      // cell_vel[0]
      //           << ", " << cell_vel[1] << ")\n";
      // vec3d u_color = Interpolators::linear(
      //     vec3d(0.0f, 0.0f, 1.0f), vec3d(1.0f, 0.0f, 0.0f), cell_vel[0] /
      //     100);

      vec3d u_color =
          cell_vel[0] > 0 ? vec3d(1.0f, 0.0f, 0.0f) : vec3d(0.0f, 0.0f, 1.0f);

      // vec3d v_color = Interpolators::linear(
      //     vec3d(0.0f, 0.0f, 1.0), vec3d(1.0f, 0.0f, 0.0), cell_vel[1] / 100);

      vec3d v_color =
          cell_vel[1] > 0 ? vec3d(1.0f, 0.0f, 0.0f) : vec3d(0.0f, 0.0f, 1.0f);

      u_color *= abs(cell_vel[0] / 10);

      v_color *= abs(cell_vel[1] / 10);

      vec3d color = 0.5f * d * (u_color + v_color);
      // color = d * v_color;
      color = {d, d, d};

      _draw_quad(x0, y0, x1, y1, color[0], color[1], color[2]);

      if (mac.get_cell_type(i, j) == CellType::SOLID) {
        // std::cout << "Cell vel of solid cell (" << i << ", " << j << "): ("
        //           << cell_vel[0] << ", " << cell_vel[1] << ")\n";
        // _draw_quad(i * mac.dx, j * mac.dy, (i + 1) * mac.dx, (j + 1) *
        // mac.dy,
        // 0.2f, 0.2f, 0.2f);
        _draw_outline(i * mac.dx, j * mac.dy, (i + 1) * mac.dx,
                      (j + 1) * mac.dy, 0.2f, 0.2f, 0.2f);
        continue;
      }
      // _draw_outline(x0, y0, x1, y1, 0.0f, 0.0f, 0.0f);
    }
  }
}

void Renderer2D::_draw_outline(float x0, float y0, float x1, float y1, float r,
                               float g, float b) {
  glBegin(GL_LINE_LOOP);
  // glColor3f(.2f, .2f, .2f);
  glColor3f(r, g, b);
  glVertex2f(x0, y0);
  glVertex2f(x1, y0);
  glVertex2f(x1, y1);
  glVertex2f(x0, y1);
  glEnd();
}

void Renderer2D::_draw_quad(float x0, float y0, float x1, float y1, float r,
                            float g, float b) {
  glColor3f(r, g, b);
  glBegin(GL_QUADS);
  glVertex2f(x0, y0);
  glVertex2f(x1, y0);
  glVertex2f(x1, y1);
  glVertex2f(x0, y1);
  glEnd();
}
