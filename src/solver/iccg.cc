#include "iccg.h"

ICCGSolver::ICCGSolver(int max_iterations, float tolerance)
    : max_iterations(max_iterations), tolerance(tolerance) {}

Eigen::VectorXd ICCGSolver::solve(const Eigen::VectorXd &divergence,
                                  const Eigen::SparseMatrix<double> &A) {
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IncompleteCholesky<double>>
      solver;

  solver.setMaxIterations(max_iterations);
  solver.setTolerance(tolerance);

  solver.compute(A);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("[SOLVER] Decomposition failed in PCG solver.");
  }

  return solver.solve(divergence);
};
