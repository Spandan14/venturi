#include <engine/mac/mac2d.h>
#include <engine/sim2d.h>
#include <renderer/renderer_2d.h>
#include <utils/physical_consts.h>

int main() {
  Simulation2D sim = Simulation2D(220, 100, 1, 1);

  auto density_init = [](int i, int j) {
    if (i >= 50 && i < 170 && j >= 50 && j < 55) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  };

  auto force_init = [](int i, int j) {
    return vec2d(0.0f, -GRAVITATIONAL_ACCL);
  };

  // auto vel_u_init = [](int i, int j) { return 3 * (j - 50); };
  //
  // auto vel_v_init = [](int i, int j) { return 0.9 * (110 - i); };

  auto vel_u_init = [](int i, int j) {
    float cx = 110, cy = 50; // center
    float dx = i - cx;
    float dy = j - cy;
    return -dy / (dx * dx + dy * dy + 1) * 500; // sharp spin
  };

  auto vel_v_init = [](int i, int j) {
    float cx = 110, cy = 50;
    float dx = i - cx;
    float dy = j - cy;
    return dx / (dx * dx + dy * dy + 1) * 500;
  };

  // auto vel_u_init = [](int i, int j) {
  //   float d1x = i - 70, d1y = j - 50;
  //   float d2x = i - 150, d2y = j - 50;
  //
  //   float f1 = -d1y / (d1x * d1x + d1y * d1y + 20) * 500;
  //   float f2 = d2y / (d2x * d2x + d2y * d2y + 20) * 500;
  //   return f1 + f2;
  // };
  //
  // auto vel_v_init = [](int i, int j) {
  //   float d1x = i - 70, d1y = j - 50;
  //   float d2x = i - 150, d2y = j - 50;
  //
  //   float f1 = d1x / (d1x * d1x + d1y * d1y + 20) * 500;
  //   float f2 = -d2x / (d2x * d2x + d2y * d2y + 20) * 500;
  //   return f1 + f2;
  // };

  sim.initialize_density(density_init);
  sim.initialize_forces(force_init);
  sim.initialize_vel_u(vel_u_init);
  sim.initialize_vel_v(vel_v_init);

  Renderer2D renderer = Renderer2D(sim.get_mac(), 1760, 800);

  auto last_frame = std::chrono::high_resolution_clock::now();
  while (renderer.should_draw()) {
    auto this_frame = std::chrono::high_resolution_clock::now();
    sim.step(std::chrono::duration<float>(this_frame - last_frame).count());
    // sim.step(0.004f);
    // std::cout << "Frame Time: "
    //           << std::chrono::duration<float, std::milli>(this_frame -
    //                                                       last_frame)
    //                  .count()
    //           << " ms" << std::endl;

    last_frame = this_frame;

    renderer.render();
    // std::cout << "Current Time: " << sim.get_mac().current_time << std::endl;
    renderer.poll();
  }
}
