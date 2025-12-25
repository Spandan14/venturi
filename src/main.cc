#include "fluxlang/runtime.h"
#include <engine/mac/mac2d.h>
#include <engine/sim2d.h>
#include <fluxlang/parser.h>
#include <fluxlang/transformer.h>
#include <peglib.h>
#include <utils/physical_consts.h>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <script.flx>" << std::endl;
    return 1;
  }

  std::string script_path = argv[1];
  Parser parser;

  // Simulation2D sim = Simulation2D(220, 100, 1, 1);
  auto ast = parser.parse(script_path.c_str());

  FluxASTTransformer transformer(*ast);
  auto script = transformer.transform();

  Runtime runtime(std::move(script));
  runtime.run();

  // auto density_init = [](int i, int j) {
  //   // if (i >= 50 && i < 170 && j >= 50 && j < 55) {
  //   if (j < 30 && j > 10 && i >= 80 && i < 140) {
  //     return 1.0f;
  //   } else {
  //     return 0.10f;
  //   }
  //   // return 0.2f;
  // };
  //
  // auto force_init = [](int i, int j) {
  //   auto force = vec2d(0.0f, -GRAVITATIONAL_ACCL);
  //
  //   float u_d1x = i - 70, u_d1y = j - 50;
  //   float u_d2x = i - 150, u_d2y = j - 50;
  //
  //   float u_f1 = -u_d1y / (u_d1x * u_d1x + u_d1y * u_d1y + 20) * 1000;
  //   float u_f2 = u_d2y / (u_d2x * u_d2x + u_d2y * u_d2y + 20) * 1000;
  //
  //   float v_d1x = i - 70, v_d1y = j - 50;
  //   float v_d2x = i - 150, v_d2y = j - 50;
  //
  //   float v_f1 = v_d1x / (v_d1x * v_d1x + v_d1y * v_d1y + 20) * 1000;
  //   float v_f2 = -v_d2x / (v_d2x * v_d2x + v_d2y * v_d2y + 20) * 1000;
  //
  //   force[0] += u_f1 + u_f2;
  //   force[1] += v_f1 + v_f2;
  //
  //   return force;
  //   // return vec2d(0.0f, 0.0f);
  // };
  //
  // auto solid_init = [](int i, int j) {
  //   // create a box in the middle
  //   if (j > 75 && j < 80 && i >= 60 && i < 160) {
  //     if (i >= 105 && i < 115) {
  //       return false;
  //     }
  //     return true;
  //   }
  //
  //   if (j > 45 && j < 50 && i >= 60 && i < 160) {
  //     if (i >= 105 && i < 115) {
  //       return false;
  //     }
  //     return true;
  //   }
  //
  //   // box around densities
  //   // if (((j < 15 && j > 8) || (j > 62 && j < 65)) && i >= 60 && i < 160) {
  //   //   return true;
  //   // }
  //   //
  //   if (((i < 60 && i > 57) || (i > 159 && i < 162)) && j > 60 && j < 70) {
  //     return true;
  //   }
  //
  //   return false;
  // };
  //
  // // auto vel_u_init = [](int i, int j) { return 0.3 * (j - 50); };
  // //
  // // auto vel_v_init = [](int i, int j) { return 0.9 * (110 - i); };
  // //
  // // auto vel_v_init = [](int i, int j) {
  // //   if (j > 70) {
  // //     return 0.5f;
  // //   } else {
  // //     return 0.0f;
  // //   }
  // // };
  //
  // // auto vel_u_init = [](int i, int j) {
  // //   float cx = 110, cy = 50; // center
  // //   float dx = i - cx;
  // //   float dy = j - cy;
  // //   return -dy / (dx * dx + dy * dy + 1) * 500; // sharp spin
  // // };
  // //
  // // auto vel_v_init = [](int i, int j) {
  // //   float cx = 110, cy = 50;
  // //   float dx = i - cx;
  // //   float dy = j - cy;
  // //   return dx / (dx * dx + dy * dy + 1) * 500;
  // // };
  //
  // auto vel_u_init = [](int i, int j) {
  //   float d1x = i - 70, d1y = j - 50;
  //   float d2x = i - 150, d2y = j - 50;
  //
  //   float f1 = -d1y / (d1x * d1x + d1y * d1y + 20) * 2000;
  //   float f2 = d2y / (d2x * d2x + d2y * d2y + 20) * 2000;
  //   return f1 + f2;
  // };
  //
  // auto vel_v_init = [](int i, int j) {
  //   float d1x = i - 70, d1y = j - 50;
  //   float d2x = i - 150, d2y = j - 50;
  //
  //   float f1 = d1x / (d1x * d1x + d1y * d1y + 20) * 1000;
  //   float f2 = -d2x / (d2x * d2x + d2y * d2y + 20) * 1000;
  //   return f1 + f2;
  // };
  //
  // sim.initialize_density(density_init);
  // sim.initialize_forces(force_init);
  // // sim.initialize_vel_u(vel_u_init);
  // // sim.initialize_vel_v(vel_v_init);
  //
  // sim.initialize_solids(solid_init);
  //
  // Renderer2D renderer = Renderer2D(sim.get_mac(), 1760, 800);
  //
  // auto last_frame = std::chrono::high_resolution_clock::now();
  // while (renderer.should_draw()) {
  //   auto this_frame = std::chrono::high_resolution_clock::now();
  //   auto step_ms =
  //       std::chrono::duration<float>(this_frame - last_frame).count();
  //   sim.step(step_ms);
  //   std::cout << "FPS: " << 1.0f / step_ms << std::endl;
  //
  //   last_frame = this_frame;
  //
  //   renderer.render();
  //
  //   renderer.poll();
  //   //
  //   // sleep for 2 seconds
  //   // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  // }
}
