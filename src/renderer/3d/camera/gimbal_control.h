#pragma once

#include <renderer/3d/camera/camera.h>

class GimbalControl {
public:
  GimbalControl(std::shared_ptr<Camera> camera, vec3f center, float radius,
                float axial_sensitivity = 0.05f,
                float scroll_sensitivity = 10.f)
      : camera(camera), center(center), radius(radius),
        axial_sensitivity(axial_sensitivity),
        scroll_sensitivity(scroll_sensitivity) {}
  ~GimbalControl() = default;

  void update_camera() {
    camera->position = vec3f{
        radius * cosf(elevation) * sinf(azimuth),
        radius * sinf(elevation),
        radius * cosf(elevation) * cosf(azimuth),
    };

    camera->position += center;

    camera->look_at(center);
  }

  void handle_mouse_move(float delta_x, float delta_y) {
    azimuth += delta_x * axial_sensitivity;
    elevation += delta_y * axial_sensitivity;

    // clamp elevation to avoid gimbal lock
    const float max_elevation = M_PI / 2.0f - 0.01f;
    if (elevation > max_elevation) {
      elevation = max_elevation;
    } else if (elevation < -max_elevation) {
      elevation = -max_elevation;
    }

    update_camera();
  }

  void handle_mouse_scroll(float delta_scroll) {
    radius -= delta_scroll * scroll_sensitivity;
    if (radius < 0.1f) {
      radius = 0.1f;
    }

    update_camera();
  }

  void reset() {
    azimuth = 0.0f;
    elevation = 0.0f;
    update_camera();
  }

private:
  std::shared_ptr<Camera> camera;
  vec3f center;
  float radius;
  float axial_sensitivity;
  float scroll_sensitivity;

  float azimuth = 0.0f;
  float elevation = 0.0f;
};
