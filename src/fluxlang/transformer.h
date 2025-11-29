#pragma once

#include "fluxlang/ast/nodes.h"
#include "peglib.h"

class FluxASTTransformer {
public:
  FluxASTTransformer(peg::Ast &root) : peg_root(root) {}
  std::unique_ptr<Flux::Script> transform();

private:
  peg::Ast &peg_root;

  std::unique_ptr<Flux::Expression> transformExpression(peg::Ast &node);
};
