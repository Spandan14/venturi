#include <fluxlang/ast/nodes.h>
#include <fluxlang/ast/visitor.h>

namespace Flux {
void Script::accept(Visitor &visitor) { visitor.visit(*this); }

void LiteralExpression::accept(Visitor &visitor) { visitor.visit(*this); }
void GenVariableExpression::accept(Visitor &visitor) { visitor.visit(*this); }
void RangeExpression::accept(Visitor &visitor) { visitor.visit(*this); }
void VectorExpression::accept(Visitor &visitor) { visitor.visit(*this); }
void CellSelectorExpression::accept(Visitor &visitor) { visitor.visit(*this); }
void BinaryExpression::accept(Visitor &visitor) { visitor.visit(*this); }
void GenFuncCallExpression::accept(Visitor &visitor) { visitor.visit(*this); }

void DimStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void GridStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void WindowStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void SetStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void DensityStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void ForceStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void SolidStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void FlowStatement::accept(Visitor &visitor) { visitor.visit(*this); }
void FlowRatioStatement::accept(Visitor &visitor) { visitor.visit(*this); }
} // namespace Flux
