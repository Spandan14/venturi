#pragma once

namespace Flux {
class Script;

class LiteralExpression;
class GenVariableExpression;
class RangeExpression;
class VectorExpression;
class CellSelectorExpression;
class BinaryExpression;
class GenFuncCallExpression;

class DimStatement;
class GridStatement;
class WindowStatement;
class SetStatement;
class DensityStatement;
class ForceStatement;
class SolidStatement;
class FlowStatement;

class Visitor {
public:
  virtual ~Visitor() = default;

  virtual void visit(Script &node) = 0;

  virtual void visit(LiteralExpression &node) = 0;
  virtual void visit(GenVariableExpression &node) = 0;
  virtual void visit(RangeExpression &node) = 0;
  virtual void visit(VectorExpression &node) = 0;
  virtual void visit(CellSelectorExpression &node) = 0;
  virtual void visit(BinaryExpression &node) = 0;
  virtual void visit(GenFuncCallExpression &node) = 0;

  virtual void visit(DimStatement &node) = 0;
  virtual void visit(GridStatement &node) = 0;
  virtual void visit(WindowStatement &node) = 0;
  virtual void visit(SetStatement &node) = 0;
  virtual void visit(DensityStatement &node) = 0;
  virtual void visit(ForceStatement &node) = 0;
  virtual void visit(SolidStatement &node) = 0;
  virtual void visit(FlowStatement &node) = 0;
};
} // namespace Flux
