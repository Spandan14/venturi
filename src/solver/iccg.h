#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

class ICCGSolver {
public:
  ICCGSolver(int max_iterations, float tolerance);
  ~ICCGSolver() = default;

  Eigen::VectorXd solve(const Eigen::VectorXd &divergence,
                        const Eigen::SparseMatrix<double> &A);

private:
  int max_iterations;
  float tolerance;
};
