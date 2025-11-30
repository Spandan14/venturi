#pragma once

#include <engine/sim.h>
#include <engine/sim2d.h>
#include <fluxlang/ast/nodes.h>
#include <fluxlang/ast/visitor.h>
#include <memory>
#include <optional>
#include <renderer/renderer.h>
#include <renderer/renderer_2d.h>

struct RuntimeConfig {
  std::optional<Flux::DimType> dim;
  std::unique_ptr<Simulation> simulation;
  std::unique_ptr<Renderer> renderer;
};

class Runtime : public Flux::Visitor {
public:
  Runtime() = default;
  ~Runtime() override = default;

  std::unique_ptr<Simulation> run();

private:
  std::unique_ptr<Flux::Script> _script;

  RuntimeConfig _config;

  void visit(Flux::Script &node) override;
  // void visit(Flux::LiteralExpression &node) override;
  // void visit(Flux::GenVariableExpression &node) override;
  // void visit(Flux::RangeExpression &node) override;
  // void visit(Flux::VectorExpression &node) override;
  // void visit(Flux::CellSelectorExpression &node) override;
  // void visit(Flux::BinaryExpression &node) override;
  // void visit(Flux::GenFuncCallExpression &node) override;
  void visit(Flux::DimStatement &node) override;
  void visit(Flux::GridStatement &node) override;
  void visit(Flux::WindowStatement &node) override;
  void visit(Flux::SetStatement &node) override;
  void visit(Flux::DensityStatement &node) override;
  void visit(Flux::ForceStatement &node) override;
  void visit(Flux::SolidStatement &node) override;
};
