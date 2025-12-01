#pragma once

#include "fluxlang/ast/nodes.h"
#include "peglib.h"
#include <memory>

class FluxASTTransformer {
public:
  FluxASTTransformer(peg::Ast &root) : peg_root(root) {}
  std::unique_ptr<Flux::Script> transform();

private:
  peg::Ast &peg_root;

  static inline bool _is_spacer(peg::Ast &node) { return node.name == "_"; }
  static peg::Ast &_reduce(peg::Ast &node);

  std::unique_ptr<Flux::Statement> _transform_statement(peg::Ast &node);
  std::unique_ptr<Flux::DimStatement> _transform_dim_statement(peg::Ast &node);
  std::unique_ptr<Flux::GridStatement>
  _transform_grid_statement(peg::Ast &node);
  std::unique_ptr<Flux::WindowStatement>
  _transform_window_statement(peg::Ast &node);
  std::unique_ptr<Flux::SetStatement> _transform_set_statement(peg::Ast &node);
  std::unique_ptr<Flux::DensityStatement>
  _transform_density_statement(peg::Ast &node);
  std::unique_ptr<Flux::ForceStatement>
  _transform_force_statement(peg::Ast &node);
  std::unique_ptr<Flux::SolidStatement>
  _transform_solid_statement(peg::Ast &node);

  std::unique_ptr<Flux::Expression> _transform_expression(peg::Ast &node);

  std::unique_ptr<Flux::BinaryExpression>
  _transform_bin_expression(peg::Ast &node, Flux::BinaryOp op);
  std::unique_ptr<Flux::BinaryExpression>
  _transform_chained_expression(peg::Ast &node);
  std::unique_ptr<Flux::BinaryExpression>
  _transform_cmp_expression(peg::Ast &node);

  std::unique_ptr<Flux::RangeExpression>
  _transform_range_expression(peg::Ast &node);

  std::unique_ptr<Flux::GenVariableExpression>
  _transform_gen_variable_expression(peg::Ast &node);
  std::unique_ptr<Flux::GenFuncCallExpression>
  _transform_gen_func_call_expression(peg::Ast &node);

  std::unique_ptr<Flux::LiteralExpression>
  _transform_literal_expression(peg::Ast &node);

  std::unique_ptr<Flux::VectorExpression>
  _transform_vector_expression(peg::Ast &node);

  std::unique_ptr<Flux::CellSelectorExpression>
  _transform_cell_selector_expression(peg::Ast &node);

  std::string _extract_target(peg::Ast &node);
  Flux::GenVar _extract_gen_var(peg::Ast &node);
  Flux::GenFunc _extract_gen_func(peg::Ast &node);
};
