#pragma once
#include <Eigen/Sparse>

namespace Preconditioners {
Eigen::VectorXf incomplete_cholesky(const Eigen::SparseMatrix<float> &A);
} // namespace Preconditioners
