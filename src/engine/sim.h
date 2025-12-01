#pragma once

#include <Eigen/Core>
#include <functional>

typedef Eigen::Vector2d vec2d;
typedef Eigen::Vector3d vec3d;

template <int DIM> class Simulation {
public:
  static_assert(DIM == 2 || DIM == 3, "Simulation dimension must be 2 or 3.");

  using DensityInitializer =
      std::conditional_t<DIM == 2, std::function<float(int, int)>,
                         std::function<float(int, int, int)>>;
  using ForceInitializer =
      std::conditional_t<DIM == 2, std::function<vec2d(int, int)>,
                         std::function<vec3d(int, int, int)>>;
  using SolidInitializer =
      std::conditional_t<DIM == 2, std::function<bool(int, int)>,
                         std::function<bool(int, int, int)>>;

  Simulation() = default;
  virtual ~Simulation() = default;

  virtual void step(float dt) = 0;

  virtual void initialize_density(const DensityInitializer &initializer) = 0;
  virtual void initialize_forces(const ForceInitializer &initializer) = 0;
  virtual void initialize_solids(const SolidInitializer &initializer) = 0;
};
