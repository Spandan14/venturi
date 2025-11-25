#include "preconditioners.h"

Eigen::VectorXf
Preconditioners::incomplete_cholesky(const Eigen::SparseMatrix<float> &A) {
  int n = A.rows();
  Eigen::VectorXf D(n);
  D.setZero();

  // Perform incomplete Cholesky factorization
  for (int i = 0; i < n; ++i) {
    float sum = 0.0f;
    for (Eigen::SparseMatrix<float>::InnerIterator it(A, i); it; ++it) {
      int j = it.col();
      if (j < i) {
        sum += A.coeff(i, j) * A.coeff(i, j) * D(j);
      }
    }
    D(i) = 1.0f / (A.coeff(i, i) - sum);
  }

  return D;
}
