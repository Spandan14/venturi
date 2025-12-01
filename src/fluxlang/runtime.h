#pragma once

#include <engine/sim.h>
#include <engine/sim2d.h>
#include <fluxlang/ast/nodes.h>
#include <fluxlang/ast/visitor.h>
#include <memory>
#include <optional>
#include <renderer/renderer.h>
#include <renderer/renderer_2d.h>

using Value2D = std::function<float(float, float)>;
using Value3D = std::function<float(float, float, float)>;
using Value = std::variant<Value2D, Value3D>;

using Membership2D = std::function<bool(int, int)>;
using Membership3D = std::function<bool(int, int, int)>;
using Membership = std::variant<Membership2D, Membership3D>;

using UnifiedValue =
    std::variant<Value, std::shared_ptr<std::vector<Value>>, Membership>;

struct RuntimeConfig {
  std::optional<Flux::DimType> dim;
  std::variant<std::monostate, std::unique_ptr<Simulation<2>>,
               std::unique_ptr<Simulation<3>>>
      simulation;
  std::unique_ptr<Renderer> renderer;

  std::map<std::string, Membership> set_memberships;
  std::vector<std::string> set_order;

  std::unordered_map<std::string, Value> density_values;
  std::unordered_map<std::string, std::shared_ptr<std::vector<Value>>>
      force_values;
  std::unordered_map<std::string, Value> solid_values;
};

class Runtime {
public:
  Runtime(std::unique_ptr<Flux::Script> script)
      : _script(std::move(script)), _config() {};
  ~Runtime() = default;

  void run();

  static std::function<float(const UnifiedValue &, float, float)>
      extract_float_2d;
  static std::function<float(const UnifiedValue &, float, float, float)>
      extract_float_3d;

private:
  std::unique_ptr<Flux::Script> _script;

  RuntimeConfig _config;

  void eval(Flux::Script &node);
  void _prepare_gen();
  void _prepare_2d();
  void _prepare_3d();
  void _step(float dt);

  UnifiedValue eval(Flux::Expression &node);
  void eval(Flux::Statement &node);

  Value eval(Flux::LiteralExpression &node);
  Value eval(Flux::GenVariableExpression &node);
  Value eval(Flux::RangeExpression &node);
  std::shared_ptr<std::vector<Value>> eval(Flux::VectorExpression &node);
  Membership eval(Flux::CellSelectorExpression &node);
  Value eval(Flux::BinaryExpression &node);
  Value eval(Flux::GenFuncCallExpression &node);

  void eval(Flux::DimStatement &node);
  void eval(Flux::GridStatement &node);
  void eval(Flux::WindowStatement &node);
  void eval(Flux::SetStatement &node);
  void eval(Flux::DensityStatement &node);
  void eval(Flux::ForceStatement &node);
  void eval(Flux::SolidStatement &node);
};
