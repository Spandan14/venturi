#include <fluxlang/ast/nodes.h>

std::string Flux::Script::to_string() {
  std::string result = "FluxLang Script {\n";
  for (const auto &stmt : statements) {
    result += "  " + stmt->to_string() + "\n";
  }
  result += "}";
  return result;
}

std::string Flux::dim_type_to_string(DimType dt) {
  switch (dt) {
  case DimType::Two:
    return "2D";
  case DimType::Three:
    return "3D";
  }
}

std::string Flux::target_type_to_string(TargetType tt) {
  switch (tt) {
  case TargetType::All:
    return "all";
  case TargetType::Identifier:
    return "identifier";
  }
}

std::string Flux::binary_op_to_string(BinaryOp op) {
  switch (op) {
  case BinaryOp::Add:
    return "+";
  case BinaryOp::Subtract:
    return "-";
  case BinaryOp::Multiply:
    return "*";
  case BinaryOp::Divide:
    return "/";
  case BinaryOp::Power:
    return "^";
  case BinaryOp::And:
    return "&&";
  case BinaryOp::Or:
    return "||";
  case BinaryOp::Equal:
    return "==";
  case BinaryOp::NotEqual:
    return "!=";
  case BinaryOp::Less:
    return "<";
  case BinaryOp::LessEqual:
    return "<=";
  case BinaryOp::Greater:
    return ">";
  case BinaryOp::GreaterEqual:
    return ">=";
  }
}

std::string Flux::gen_var_to_string(GenVar gv) {
  switch (gv) {
  case GenVar::I:
    return "i";
  case GenVar::J:
    return "j";
  case GenVar::K:
    return "k";
  case GenVar::T:
    return "t";
  }
}

std::string Flux::gen_func_to_string(GenFunc gf) {
  switch (gf) {
  case GenFunc::Sin:
    return "sin";
  case GenFunc::Cos:
    return "cos";
  case GenFunc::Tan:
    return "tan";
  case GenFunc::Abs:
    return "abs";
  case GenFunc::Sqrt:
    return "sqrt";
  case GenFunc::Log:
    return "log";
  case GenFunc::Exp:
    return "exp";
  }
}

// EXPRESSIONS
std::string Flux::LiteralExpression::to_string() {
  return std::to_string(value);
}

std::string Flux::GenVariableExpression::to_string() {
  return gen_var_to_string(var);
}

std::string Flux::RangeExpression::to_string() {
  return "range(" + start->to_string() + " <= " + gen_var_to_string(var) +
         " <= " + end->to_string() + ")";
}

std::string Flux::VectorExpression::to_string() {
  std::string result = "vec(";
  for (size_t i = 0; i < components.size(); ++i) {
    result += components[i]->to_string();
    if (i < components.size() - 1) {
      result += ", ";
    }
  }
  result += ")";
  return result;
}

std::string Flux::CellSelectorExpression::to_string() {
  return "cell_selector(" + condition->to_string() + ")";
}

std::string Flux::BinaryExpression::to_string() {
  return "(" + left->to_string() + " " + binary_op_to_string(op) + " " +
         right->to_string() + ")";
}

std::string Flux::GenFuncCallExpression::to_string() {
  return gen_func_to_string(func) + "(" + argument->to_string() + ")";
}

// STATEMENTS
std::string Flux::DimStatement::to_string() {
  return "dim " + dim_type_to_string(dim) + ";";
}

std::string Flux::GridStatement::to_string() {
  std::string result = "grid (";
  for (size_t i = 0; i < sizes.size(); ++i) {
    result += sizes[i]->to_string();
    if (i < sizes.size() - 1) {
      result += ", ";
    }
  }
  result += ") with dx ";
  result += std::to_string(dx) + ";";
  return result;
}

std::string Flux::WindowStatement::to_string() {
  std::string result = "window (";
  for (size_t i = 0; i < sizes.size(); ++i) {
    result += sizes[i]->to_string();
    if (i < sizes.size() - 1) {
      result += ", ";
    }
  }
  result += ");";
  return result;
}

std::string Flux::SetStatement::to_string() {
  return "set " + identifier + " " + selector->to_string() + ";";
}

std::string Flux::DensityStatement::to_string() {
  return "density " + identifier + " = " + value->to_string() + ";";
}

std::string Flux::ForceStatement::to_string() {
  return "force " + identifier + " = " + value->to_string() + ";";
}

std::string Flux::SolidStatement::to_string() {
  return "solid " + identifier + " = " + value->to_string() + ";";
}

std::string Flux::FlowStatement::to_string() {
  return "flow " + identifier + " = " + value->to_string() + ";";
}
