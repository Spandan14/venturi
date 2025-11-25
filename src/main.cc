#include <engine/mac/mac2d.h>
#include <engine/sim2d.h>
#include <iostream>
#include <renderer/renderer_2d.h>

int main() {
  Simulation2D sim = Simulation2D(220, 100, 0.1, 0.1);

  auto density_init = [](int i, int j) {
    if (i >= 55 && i < 80 && j >= 50 && j < 55) {
      return 1.0f;
    } else {
      return 0.0f;
    }
  };

  auto force_init = [](int i, int j) { return vec2d(0.0f, -9.81f); };

  auto vel_u_init = [](int i, int j) { return 0.3 * (j - 50); };

  auto vel_v_init = [](int i, int j) { return 0.15 * (110 - i); };

  sim.initialize_density(density_init);
  sim.initialize_forces(force_init);
  sim.initialize_vel_u(vel_u_init);
  sim.initialize_vel_v(vel_v_init);

  Renderer2D renderer = Renderer2D(sim.get_mac(), 1760, 800);

  auto last_frame = std::chrono::high_resolution_clock::now();
  while (renderer.should_draw()) {
    auto this_frame = std::chrono::high_resolution_clock::now();
    sim.step(std::chrono::duration<float>(this_frame - last_frame).count());
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
