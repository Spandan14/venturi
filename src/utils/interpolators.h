#pragma once

namespace Interpolators {

// requires x, y in [0, 1]
template <typename T> T bilinear(T q11, T q12, T q21, T q22, float x, float y) {
  return (q11 * (1 - x) * (1 - y) + q21 * x * (1 - y) + q12 * (1 - x) * y +
          q22 * x * y);
}

// requires x, y, z in [0, 1]
template <typename T>
T trilinear(T q000, T q001, T q010, T q011, T q100, T q101, T q110, T q111,
            float x, float y, float z) {
  return (q000 * (1 - x) * (1 - y) * (1 - z) + q100 * x * (1 - y) * (1 - z) +
          q010 * (1 - x) * y * (1 - z) + q110 * x * y * (1 - z) +
          q001 * (1 - x) * (1 - y) * z + q101 * x * (1 - y) * z +
          q011 * (1 - x) * y * z + q111 * x * y * z);
}

} // namespace Interpolators
