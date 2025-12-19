#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

typedef Eigen::Vector3f vec3f;
typedef Eigen::Quaternionf quatf;
typedef Eigen::Matrix3f mat3f;
typedef Eigen::Matrix4f mat4f;

const float DEFAULT_FOV = 60.0;
const float DEFAULT_ASPECT_RATIO = 1.0;
const float DEFAULT_NEAR_CLIP = 0.1;
const float DEFAULT_FAR_CLIP = 1000.0;

class Camera {
public:
  Camera(vec3f pos = vec3f(0.0, 0.0, 0.0), quatf orient = quatf::Identity(),
         float fov_deg = DEFAULT_FOV, float aspect = DEFAULT_ASPECT_RATIO,
         float near_c = DEFAULT_NEAR_CLIP, float far_c = DEFAULT_FAR_CLIP)
      : position(pos), orientation(orient), fov(fov_deg), aspect_ratio(aspect),
        near_clip(near_c), far_clip(far_c) {}
  ~Camera() = default;

  vec3f position;
  quatf orientation;

  float fov = DEFAULT_FOV;
  float aspect_ratio = DEFAULT_ASPECT_RATIO;
  float near_clip = DEFAULT_NEAR_CLIP;
  float far_clip = DEFAULT_FAR_CLIP;

  vec3f forward() const { return orientation * vec3f(0.0, 0.0, -1.0); }

  vec3f up() const { return orientation * vec3f(0.0, 1.0, 0.0); }

  vec3f right() const { return orientation * vec3f(1.0, 0.0, 0.0); }

  mat4f get_view_matrix() const {
    mat3f rot_matrix = orientation.conjugate().toRotationMatrix();
    mat4f view = mat4f::Identity();

    view.block<3, 3>(0, 0) = rot_matrix;
    view.block<3, 1>(0, 3) = -rot_matrix * position;

    return view;
  }

  mat4f get_view_matrix_inverse() const {
    mat3f rot_matrix = orientation.toRotationMatrix();
    mat4f inv_view = mat4f::Identity();

    inv_view.block<3, 3>(0, 0) = rot_matrix;
    inv_view.block<3, 1>(0, 3) = position;

    return inv_view;
  }

  mat4f get_projection_matrix() const {
    float fov_rad = fov * M_PI / 180.0;
    float f = 1.0 / tan(fov_rad / 2.0);

    mat4f proj = mat4f::Zero();
    proj(0, 0) = f / aspect_ratio;
    proj(1, 1) = f;
    proj(2, 2) = (far_clip + near_clip) / (near_clip - far_clip);
    proj(2, 3) = (2.0 * far_clip * near_clip) / (near_clip - far_clip);
    proj(3, 2) = -1.0;

    return proj;
  }

  mat4f get_projection_matrix_inverse() const {
    float fov_rad = fov * M_PI / 180.0;
    float f = 1.0 / tan(fov_rad / 2.0);

    mat4f inv_proj = mat4f::Zero();
    inv_proj(0, 0) = aspect_ratio / f;
    inv_proj(1, 1) = 1.0 / f;
    inv_proj(2, 3) = (near_clip - far_clip) / (2.0 * far_clip * near_clip);
    inv_proj(3, 2) = -1.0;
    inv_proj(3, 3) = (far_clip + near_clip) / (2.0 * far_clip * near_clip);

    return inv_proj;
  }
};
