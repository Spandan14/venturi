#include "fluxlang/ast/nodes.h"
#include <fluxlang/runtime.h>

std::unique_ptr<Simulation> Runtime::run() {
  _script->accept(*this);
  return std::move(_config.simulation);
}

void Runtime::visit(Flux::Script &node) {
  for (auto &stmt : node.statements) {
    stmt->accept(*this);
  }
};

void Runtime::visit(Flux::DimStatement &node) {
  _config.dim = node.dim;
  if (_config.dim != Flux::DimType::Two ||
      _config.dim != Flux::DimType::Three) {
    throw std::runtime_error("Unsupported dimension type!");
  }
};

void Runtime::visit(Flux::GridStatement &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before GridStatement!");
  }

  if (_config.dim == Flux::DimType::Two) {
    if (node.sizes.size() != 2) {
      throw std::runtime_error(
          "GridStatement size mismatch for 2D simulation!");
    }

    int nx = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[0].get())->value);
    int ny = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[1].get())->value);

    _config.simulation = std::make_unique<Simulation2D>(
        nx, ny, static_cast<float>(node.dx), static_cast<float>(node.dx));
  }

  if (_config.dim == Flux::DimType::Three) {
    throw std::runtime_error("3D simulations are not yet supported!");
  }
};

void Runtime::visit(Flux::WindowStatement &node) {
  if (!_config.dim.has_value()) {
    throw std::runtime_error(
        "Dimension type must be specified before WindowStatement!");
  }

  if (_config.dim == Flux::DimType::Two) {
    if (node.sizes.size() != 2) {
      throw std::runtime_error(
          "WindowStatement size mismatch for 2D simulation!");
    }

    int width = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[0].get())->value);
    int height = static_cast<int>(
        dynamic_cast<Flux::LiteralExpression *>(node.sizes[1].get())->value);

    // Create renderer
    _config.renderer = std::make_unique<Renderer2D>(
        dynamic_cast<Simulation2D *>(_config.simulation.get())->get_mac(),
        width, height);
  }

  if (_config.dim == Flux::DimType::Three) {
    throw std::runtime_error("3D simulations are not yet supported!");
  }
};

void Runtime::visit(Flux::SetStatement &node) {};

void Runtime::visit(Flux::DensityStatement &node) {};

void Runtime::visit(Flux::ForceStatement &node) {};

void Runtime::visit(Flux::SolidStatement &node) {};
