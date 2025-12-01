#include <fluxlang/transformer.h>

std::unique_ptr<Flux::Script> FluxASTTransformer::transform() {
  // first child of the root is always a Script
  if (peg_root.name != "Script") {
    throw std::runtime_error("Expected Script node as first child of root!");
  }

  std::vector<std::unique_ptr<Flux::Statement>> statements;

  for (auto &stmt_child : peg_root.nodes) {
    if (_is_spacer(*stmt_child)) {
      continue;
    }

    auto stmt = _transform_statement(*stmt_child);
    std::cout << "Transformed Statement: " << stmt->to_string() << std::endl;
    statements.push_back(std::move(stmt));
  }

  return std::make_unique<Flux::Script>(std::move(statements));
}

peg::Ast &FluxASTTransformer::_reduce(peg::Ast &node) {
  // recursively reduce child nodes
  if (node.nodes.size() != 1) {
    return node;
  }

  return _reduce(*node.nodes[0]);
}

std::unique_ptr<Flux::Statement>
FluxASTTransformer::_transform_statement(peg::Ast &node) {
  if (node.nodes.empty()) {
    throw std::runtime_error("Statement node has no children!");
  }

  if (node.nodes.size() != 1) {
    throw std::runtime_error("Statement node has more than one child!");
  }

  auto &stmt_node = node.nodes[0];

  // FIXME: clean this up somehow ?
  if (stmt_node->name == "DimStmt") {
    return _transform_dim_statement(*stmt_node);
  } else if (stmt_node->name == "GridStmt") {
    return _transform_grid_statement(*stmt_node);
  } else if (stmt_node->name == "WindowStmt") {
    return _transform_window_statement(*stmt_node);
  } else if (stmt_node->name == "SetStmt") {
    return _transform_set_statement(*stmt_node);
  } else if (stmt_node->name == "DensityStmt") {
    return _transform_density_statement(*stmt_node);
  } else if (stmt_node->name == "ForceStmt") {
    return _transform_force_statement(*stmt_node);
  } else if (stmt_node->name == "SolidStmt") {
    return _transform_solid_statement(*stmt_node);
  } else {
    throw std::runtime_error("Unknown Statement type: " + stmt_node->name);
  }
};

std::unique_ptr<Flux::DimStatement>
FluxASTTransformer::_transform_dim_statement(peg::Ast &node) {
  // DimStmt        <- 'dim' _ DimType
  // look for a DimType, if not found, then throw parsing error

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "DimType") {
      if (!child->is_token) {
        throw std::runtime_error("DimType node is not a token!");
      }

      std::string token = std::string(child->token);

      if (token == "2D") {
        return std::make_unique<Flux::DimStatement>(Flux::DimType::Two);
      } else if (token == "3D") {
        return std::make_unique<Flux::DimStatement>(Flux::DimType::Three);
      } else {
        throw std::runtime_error("Unknown DimType token: " + token);
      }
    }
  }

  // unreachable
  throw std::runtime_error("DimType not found in DimStmt!");
}

std::unique_ptr<Flux::GridStatement>
FluxASTTransformer::_transform_grid_statement(peg::Ast &node) {
  // GridStmt       <- 'grid' _ Int _ 'x' _ Int ( _ 'x' _ Int )?
  std::vector<int> sizes;
  std::optional<double> dx = std::nullopt;
  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Int") {
      if (!child->is_token) {
        throw std::runtime_error("Int node is not a token!");
      }

      int size = std::stoi(std::string(child->token));
      sizes.push_back(size);
    }

    if (child->name == "Number") {
      if (!child->is_token) {
        throw std::runtime_error("Number node is not a token!");
      }

      dx = std::stod(std::string(child->token));
    }
  }

  if (sizes.size() != 2 && sizes.size() != 3) {
    throw std::runtime_error("GridStmt must have 2 or 3 sizes!");
  }

  if (!dx.has_value()) {
    throw std::runtime_error("GridStmt must have a dx value!");
  }

  std::vector<std::unique_ptr<Flux::Expression>> size_exprs;
  for (int size : sizes) {
    size_exprs.push_back(
        std::make_unique<Flux::LiteralExpression>(static_cast<double>(size)));
  }

  return std::make_unique<Flux::GridStatement>(std::move(size_exprs),
                                               dx.value());
}

std::unique_ptr<Flux::WindowStatement>
FluxASTTransformer::_transform_window_statement(peg::Ast &node) {
  // WindowStmt     <- 'window' _ Int _ 'x' _ Int ( _ 'x' _ Int )?
  std::vector<int> sizes;
  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Int") {
      if (!child->is_token) {
        throw std::runtime_error("Int node is not a token!");
      }

      int size = std::stoi(std::string(child->token));
      sizes.push_back(size);
    }
  }

  if (sizes.size() != 2 && sizes.size() != 3) {
    throw std::runtime_error("WindowStmt must have 2 or 3 sizes!");
  }

  std::vector<std::unique_ptr<Flux::Expression>> size_exprs;
  for (int size : sizes) {
    size_exprs.push_back(
        std::make_unique<Flux::LiteralExpression>(static_cast<double>(size)));
  }

  return std::make_unique<Flux::WindowStatement>(std::move(size_exprs));
}

std::unique_ptr<Flux::SetStatement>
FluxASTTransformer::_transform_set_statement(peg::Ast &node) {
  // SetStmt        <- 'set' _ Identifier _ CellSelector
  std::string identifier;
  std::unique_ptr<Flux::CellSelectorExpression> selector;
  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Identifier") {
      if (!child->is_token) {
        throw std::runtime_error("Identifier node is not a token!");
      }

      identifier = std::string(child->token);
    } else if (child->name == "CellSelector") {
      selector = _transform_cell_selector_expression(*child);
    }
  }

  if (identifier.empty()) {
    throw std::runtime_error("Identifier not found in SetStmt!");
  }

  if (!selector) {
    throw std::runtime_error("CellSelector not found in SetStmt!");
  }

  return std::make_unique<Flux::SetStatement>(std::move(identifier),
                                              std::move(selector));
};

std::unique_ptr<Flux::DensityStatement>
FluxASTTransformer::_transform_density_statement(peg::Ast &node) {
  // DensityStmt   <- 'density' _ TargetType _ Identifier _ '=' _ Expression
  Flux::TargetType target_type = Flux::TargetType::All;
  std::string identifier;
  std::unique_ptr<Flux::Expression> value;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Target") {
      identifier = _extract_target(*child);
      if (identifier != "all") {
        target_type = Flux::TargetType::Identifier;
      }
    } else if (child->name == "Expr") {
      value = _transform_expression(*child);
    }
  }

  if (identifier.empty()) {
    throw std::runtime_error("Identifier not found in DensityStmt!");
  }

  if (!value) {
    throw std::runtime_error("Expression not found in DensityStmt!");
  }

  return std::make_unique<Flux::DensityStatement>(
      target_type, std::move(identifier), std::move(value));
}

std::unique_ptr<Flux::ForceStatement>
FluxASTTransformer::_transform_force_statement(peg::Ast &node) {
  // ForceStmt     <- 'force' _ TargetType _ Identifier _ '=' _ VectorExpr
  Flux::TargetType target_type = Flux::TargetType::All;
  std::string identifier;
  std::unique_ptr<Flux::VectorExpression> value;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Target") {
      identifier = _extract_target(*child);
      if (identifier != "all") {
        target_type = Flux::TargetType::Identifier;
      }
    } else if (child->name == "VectorExpr") {
      value = _transform_vector_expression(*child);
    }
  }

  if (identifier.empty()) {
    throw std::runtime_error("Identifier not found in ForceStmt!");
  }
  if (!value) {
    throw std::runtime_error("VectorExpression not found in ForceStmt!");
  }
  return std::make_unique<Flux::ForceStatement>(
      target_type, std::move(identifier), std::move(value));
}

std::unique_ptr<Flux::SolidStatement>
FluxASTTransformer::_transform_solid_statement(peg::Ast &node) {
  // SolidStmt     <- 'solid' _ TargetType _ Identifier _ '=' _ Expression
  Flux::TargetType target_type = Flux::TargetType::All;
  std::string identifier;
  std::unique_ptr<Flux::Expression> value;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Target") {
      identifier = _extract_target(*child);
      if (identifier != "all") {
        target_type = Flux::TargetType::Identifier;
      }
    } else if (child->name == "Expr") {
      value = _transform_expression(*child);
    }
  }

  if (identifier.empty()) {
    throw std::runtime_error("Identifier not found in SolidStmt!");
  }
  if (!value) {
    throw std::runtime_error("Expression not found in SolidStmt!");
  }
  return std::make_unique<Flux::SolidStatement>(
      target_type, std::move(identifier), std::move(value));
}

std::unique_ptr<Flux::Expression>
FluxASTTransformer::_transform_expression(peg::Ast &node) {
  // Expr          <- OrExpr
  auto &reduced_expr = _reduce(node);

  if (reduced_expr.name == "OrExpr") {
    return _transform_bin_expression(reduced_expr, Flux::BinaryOp::Or);
  } else if (reduced_expr.name == "AndExpr") {
    return _transform_bin_expression(reduced_expr, Flux::BinaryOp::And);
  } else if (reduced_expr.name == "CmpExpr") {
    return _transform_cmp_expression(reduced_expr);
  } else if (reduced_expr.name == "AddExpr") {
    return _transform_chained_expression(reduced_expr);
  } else if (reduced_expr.name == "MulExpr") {
    return _transform_chained_expression(reduced_expr);
  } else if (reduced_expr.name == "Number" || reduced_expr.name == "Boolean") {
    return _transform_literal_expression(reduced_expr);
  } else if (reduced_expr.name == "Var") {
    return _transform_gen_variable_expression(reduced_expr);
  } else if (reduced_expr.name == "RangeExpr") {
    return _transform_range_expression(reduced_expr);
  } else if (reduced_expr.name == "FuncCall") {
    return _transform_gen_func_call_expression(reduced_expr);
  } else if (reduced_expr.name == "VectorExpr") {
    return _transform_vector_expression(reduced_expr);
  } else if (reduced_expr.name == "CellSelector") {
    return _transform_cell_selector_expression(reduced_expr);
  } else if (reduced_expr.name == "Primary") {
    for (auto &child : reduced_expr.nodes) {
      if (_is_spacer(*child)) {
        continue;
      }
      return _transform_expression(*child);
    }
    throw std::runtime_error("Primary node has no valid children!");
  } else {
    throw std::runtime_error("Unknown Expression type: " + reduced_expr.name);
  }
}

std::unique_ptr<Flux::BinaryExpression>
FluxASTTransformer::_transform_bin_expression(peg::Ast &node,
                                              Flux::BinaryOp op) {
  std::vector<std::unique_ptr<Flux::Expression>> operands;
  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    auto expr = _transform_expression(*child);
    operands.push_back(std::move(expr));
  }

  if (operands.size() != 2) {
    throw std::runtime_error(
        "BinaryExpression must have exactly two operands!");
  }

  return std::make_unique<Flux::BinaryExpression>(op, std::move(operands[0]),
                                                  std::move(operands[1]));
}

std::unique_ptr<Flux::BinaryExpression>
FluxASTTransformer::_transform_chained_expression(peg::Ast &node) {
  std::vector<std::unique_ptr<Flux::Expression>> operands;
  std::vector<Flux::BinaryOp> ops;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "MulOp" || child->name == "AddOp") {
      if (!child->is_token) {
        throw std::runtime_error("Operator node is not a token!");
      }

      std::string token = std::string(child->token);
      Flux::BinaryOp op;
      if (token == "+") {
        op = Flux::BinaryOp::Add;
      } else if (token == "-") {
        op = Flux::BinaryOp::Subtract;
      } else if (token == "*") {
        op = Flux::BinaryOp::Multiply;
      } else if (token == "/") {
        op = Flux::BinaryOp::Divide;
      } else if (token == "^") {
        op = Flux::BinaryOp::Power;
      } else {
        throw std::runtime_error("Unknown operator: " + token);
      }
      ops.push_back(op);
    } else {
      auto expr = _transform_expression(*child);
      operands.push_back(std::move(expr));
    }
  }

  if (operands.size() < 2) {
    throw std::runtime_error(
        "ChainedExpression must have at least two operands!");
  }

  if (ops.size() != operands.size() - 1) {
    std::cout << "Operands size: " << operands.size()
              << ", Ops size: " << ops.size() << std::endl;
    for (const auto &op : ops) {
      std::cout << "Op: " << static_cast<int>(op) << std::endl;
    }
    for (const auto &operand : operands) {
      std::cout << "Operand: " << operand->to_string() << std::endl;
    }
    throw std::runtime_error(
        "Number of operators must be one less than number of operands!");
  }

  // build left-associative binary expressions
  std::unique_ptr<Flux::Expression> left = std::move(operands[0]);
  for (size_t i = 0; i < ops.size(); ++i) {
    left = std::make_unique<Flux::BinaryExpression>(ops[i], std::move(left),
                                                    std::move(operands[i + 1]));
  }

  return std::unique_ptr<Flux::BinaryExpression>(
      static_cast<Flux::BinaryExpression *>(left.release()));
}

std::unique_ptr<Flux::BinaryExpression>
FluxASTTransformer::_transform_cmp_expression(peg::Ast &node) {
  std::vector<std::unique_ptr<Flux::Expression>> operands;
  Flux::BinaryOp op;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "CmpOp") {
      if (!child->is_token) {
        throw std::runtime_error("CmpOp node is not a token!");
      }

      std::string token = std::string(child->token);
      if (token == "==") {
        op = Flux::BinaryOp::Equal;
      } else if (token == "!=") {
        op = Flux::BinaryOp::NotEqual;
      } else if (token == "<") {
        op = Flux::BinaryOp::Less;
      } else if (token == "<=") {
        op = Flux::BinaryOp::LessEqual;
      } else if (token == ">") {
        op = Flux::BinaryOp::Greater;
      } else if (token == ">=") {
        op = Flux::BinaryOp::GreaterEqual;
      } else {
        throw std::runtime_error("Unknown comparison operator: " + token);
      }
    } else {
      auto expr = _transform_expression(*child);
      operands.push_back(std::move(expr));
    }
  }

  if (operands.size() != 2) {
    throw std::runtime_error(
        "BinaryExpression must have exactly two operands!");
  }

  return std::make_unique<Flux::BinaryExpression>(op, std::move(operands[0]),
                                                  std::move(operands[1]));
}

std::unique_ptr<Flux::RangeExpression>
FluxASTTransformer::_transform_range_expression(peg::Ast &node) {
  // RangeExpr     <- Expr _ GenVar _ Expr
  std::unique_ptr<Flux::Expression> start;
  Flux::GenVar var;
  std::unique_ptr<Flux::Expression> end;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "AddExpr") {
      if (!start) {
        start = _transform_expression(*child);
      } else {
        end = _transform_expression(*child);
      }
    } else if (child->name == "Var") {
      var = _extract_gen_var(*child);
    }
  }

  if (!start) {
    throw std::runtime_error("Start expression not found in RangeExpr!");
  }
  if (!end) {
    throw std::runtime_error("End expression not found in RangeExpr!");
  }

  return std::make_unique<Flux::RangeExpression>(std::move(start), var,
                                                 std::move(end));
}

std::unique_ptr<Flux::GenVariableExpression>
FluxASTTransformer::_transform_gen_variable_expression(peg::Ast &node) {
  // Var           <- 'i' / 'j' / 'k'
  Flux::GenVar var = _extract_gen_var(node);
  return std::make_unique<Flux::GenVariableExpression>(var);
}

std::unique_ptr<Flux::GenFuncCallExpression>
FluxASTTransformer::_transform_gen_func_call_expression(peg::Ast &node) {
  // FuncCall      <- GenFunc _ '(' _ Expr _ ')'
  Flux::GenFunc func;
  std::unique_ptr<Flux::Expression> argument;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Identifier") {
      func = _extract_gen_func(*child);
    } else if (child->name == "Expr") {
      argument = _transform_expression(*child);
    }
  }

  if (!argument) {
    throw std::runtime_error("Argument expression not found in FuncCall!");
  }

  return std::make_unique<Flux::GenFuncCallExpression>(func,
                                                       std::move(argument));
}

std::unique_ptr<Flux::LiteralExpression>
FluxASTTransformer::_transform_literal_expression(peg::Ast &node) {
  // Number        <- [0-9]+ ('.' [0-9]+)?
  // Boolean       <- 'true' / 'false'
  if (node.name == "Number") {
    if (!node.is_token) {
      throw std::runtime_error("Number node is not a token!");
    }

    double value = std::stod(std::string(node.token));
    return std::make_unique<Flux::LiteralExpression>(value);
  } else if (node.name == "Boolean") {
    if (!node.is_token) {
      throw std::runtime_error("Boolean node is not a token!");
    }

    std::string token = std::string(node.token);
    double value = (token == "true") ? 1.0 : 0.0;
    return std::make_unique<Flux::LiteralExpression>(value);
  } else {
    throw std::runtime_error("Unknown LiteralExpression type: " + node.name);
  }
}

std::unique_ptr<Flux::VectorExpression>
FluxASTTransformer::_transform_vector_expression(peg::Ast &node) {
  // VectorExpr    <- '[' _ Expr ( _ ',' _ Expr )* _ ']'
  std::vector<std::unique_ptr<Flux::Expression>> components;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Expr") {
      auto expr = _transform_expression(*child);
      components.push_back(std::move(expr));
    }
  }

  if (components.empty()) {
    throw std::runtime_error(
        "VectorExpression must have at least one component!");
  }

  return std::make_unique<Flux::VectorExpression>(std::move(components));
}

std::unique_ptr<Flux::CellSelectorExpression>
FluxASTTransformer::_transform_cell_selector_expression(peg::Ast &node) {
  // CellSelector <- 'cells' _ 'where' _ Expr
  std::unique_ptr<Flux::Expression> condition;

  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Expr") {
      condition = _transform_expression(*child);
    }
  }

  if (!condition) {
    throw std::runtime_error("Condition expression not found in CellSelector!");
  }

  return std::make_unique<Flux::CellSelectorExpression>(std::move(condition));
}

std::string FluxASTTransformer::_extract_target(peg::Ast &node) {
  std::string identifier = "all";
  for (auto &child : node.nodes) {
    if (_is_spacer(*child)) {
      continue;
    }

    if (child->name == "Identifier") {
      if (!child->is_token) {
        throw std::runtime_error("Identifier node is not a token!");
      }

      identifier = std::string(child->token);
      return identifier;
    }
  }

  return identifier;
}

Flux::GenVar FluxASTTransformer::_extract_gen_var(peg::Ast &node) {
  if (!node.is_token) {
    throw std::runtime_error("GenVar node is not a token!");
  }

  std::string token = std::string(node.token);
  if (token == "i") {
    return Flux::GenVar::I;
  } else if (token == "j") {
    return Flux::GenVar::J;
  } else if (token == "k") {
    return Flux::GenVar::K;
  } else {
    throw std::runtime_error("Unknown GenVar token: " + token);
  }
}

Flux::GenFunc FluxASTTransformer::_extract_gen_func(peg::Ast &node) {
  if (!node.is_token) {
    throw std::runtime_error("GenFunc node is not a token!");
  }

  std::string token = std::string(node.token);
  if (token == "sin") {
    return Flux::GenFunc::Sin;
  } else if (token == "cos") {
    return Flux::GenFunc::Cos;
  } else if (token == "tan") {
    return Flux::GenFunc::Tan;
  } else if (token == "abs") {
    return Flux::GenFunc::Abs;
  } else if (token == "sqrt") {
    return Flux::GenFunc::Sqrt;
  } else if (token == "log") {
    return Flux::GenFunc::Log;
  } else if (token == "exp") {
    return Flux::GenFunc::Exp;
  } else {
    throw std::runtime_error("Unknown GenFunc token: " + token);
  }
}
