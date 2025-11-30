#pragma once

#include <memory>
#include <string>
#include <vector>

namespace Flux {
class Visitor;

class ASTNode {
public:
  virtual ~ASTNode() = default;
  virtual void accept(Visitor &visitor) = 0;
  virtual std::string to_string() = 0;
};

class Expression : public ASTNode {
public:
  virtual ~Expression() = default;
};
class Statement : public ASTNode {
public:
  virtual ~Statement() = default;
};

// SCRIPT
class Script : public ASTNode {
public:
  std::vector<std::unique_ptr<Statement>> statements;
  Script(std::vector<std::unique_ptr<Statement>> stmts)
      : statements(std::move(stmts)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

// ENUMS
enum class DimType { Two, Three };

enum class TargetType { All, Identifier };

enum class BinaryOp {
  Add,
  Subtract,
  Multiply,
  Divide,
  And,
  Or,
  Equal,
  NotEqual,
  Less,
  LessEqual,
  Greater,
  GreaterEqual
};

enum class GenVar { I, J, K };

enum class GenFunc { Sin, Cos, Tan, Abs, Sqrt, Log, Exp };

std::string dim_type_to_string(DimType dt);
std::string target_type_to_string(TargetType tt);
std::string binary_op_to_string(BinaryOp op);
std::string gen_var_to_string(GenVar gv);
std::string gen_func_to_string(GenFunc gf);

// EXPRESSIONS
class LiteralExpression : public Expression {
public:
  double value;
  LiteralExpression(double val) : value(val) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class GenVariableExpression : public Expression {
public:
  GenVar var;
  GenVariableExpression(GenVar v) : var(v) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class RangeExpression : public Expression {
public:
  std::unique_ptr<Expression> start;
  GenVar var;
  std::unique_ptr<Expression> end;
  RangeExpression(std::unique_ptr<Expression> s, GenVar v,
                  std::unique_ptr<Expression> e)
      : start(std::move(s)), var(v), end(std::move(e)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class VectorExpression : public Expression {
public:
  std::vector<std::unique_ptr<Expression>> components;
  VectorExpression(std::vector<std::unique_ptr<Expression>> a)
      : components(std::move(a)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class CellSelectorExpression : public Expression {
public:
  std::unique_ptr<Expression> condition;
  CellSelectorExpression(std::unique_ptr<Expression> cond)
      : condition(std::move(cond)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class BinaryExpression : public Expression {
public:
  BinaryOp op;
  std::unique_ptr<Expression> left;
  std::unique_ptr<Expression> right;
  BinaryExpression(BinaryOp oper, std::unique_ptr<Expression> lhs,
                   std::unique_ptr<Expression> rhs)
      : op(oper), left(std::move(lhs)), right(std::move(rhs)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class GenFuncCallExpression : public Expression {
public:
  GenFunc func;
  std::unique_ptr<Expression> argument;
  GenFuncCallExpression(GenFunc f, std::unique_ptr<Expression> arg)
      : func(f), argument(std::move(arg)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

// STATEMENTS
class DimStatement : public Statement {
public:
  DimType dim;
  DimStatement(DimType d) : dim(d) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class GridStatement : public Statement {
public:
  std::vector<std::unique_ptr<Expression>> sizes;
  GridStatement(std::vector<std::unique_ptr<Expression>> r)
      : sizes(std::move(r)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class WindowStatement : public Statement {
public:
  std::vector<std::unique_ptr<Expression>> sizes;
  WindowStatement(std::vector<std::unique_ptr<Expression>> r)
      : sizes(std::move(r)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class SetStatement : public Statement {
public:
  std::string identifier;
  std::unique_ptr<CellSelectorExpression> selector;
  SetStatement(std::string id, std::unique_ptr<CellSelectorExpression> sel)
      : identifier(std::move(id)), selector(std::move(sel)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class DensityStatement : public Statement {
public:
  TargetType target_type;
  std::string identifier;
  std::unique_ptr<Expression> value;
  DensityStatement(TargetType tt, std::string id,
                   std::unique_ptr<Expression> val)
      : target_type(tt), identifier(std::move(id)), value(std::move(val)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class ForceStatement : public Statement {
public:
  TargetType target_type;
  std::string identifier;
  std::unique_ptr<VectorExpression> value;
  ForceStatement(TargetType tt, std::string id,
                 std::unique_ptr<VectorExpression> val)
      : target_type(tt), identifier(std::move(id)), value(std::move(val)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

class SolidStatement : public Statement {
public:
  TargetType target_type;
  std::string identifier;
  std::unique_ptr<Expression> value;
  SolidStatement(TargetType tt, std::string id, std::unique_ptr<Expression> val)
      : target_type(tt), identifier(std::move(id)), value(std::move(val)) {}
  void accept(Visitor &visitor) override;
  std::string to_string() override;
};

} // namespace Flux
