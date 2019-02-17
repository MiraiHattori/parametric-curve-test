#include <iostream>
#include <OsqpEigen/OsqpEigen.h>
#include <Eigen/Dense>
#include <chrono>

/*
 * minimize 1/2 x^T H x + x^T g
 * subject to lb <= x <= ub
 *            lbA <= Ax <= ubA
 */

namespace hrp
{
    using dmatrix = Eigen::MatrixXd;
    using dvector = Eigen::VectorXd;
}

// {{{ osqpWrapper
double solve_qp(
        // const hrp::dvector& initial_state,
        const hrp::dvector& state_min_vector,
        const hrp::dvector& state_max_vector,
        const hrp::dmatrix& eval_weight_matrix,
        const hrp::dvector& eval_coeff_vector,
        const hrp::dmatrix& equality_matrix,
        const hrp::dvector& equality_coeff_vector,
        const hrp::dmatrix& inequality_matrix,
        const hrp::dvector& inequality_min_vector,
        const hrp::dvector& inequality_max_vector,
        hrp::dvector& result_vector
        )
{
  int x_size = static_cast<int>(eval_weight_matrix.rows());
  int eq_size = static_cast<int>(equality_coeff_vector.size());
  int ieq_size = static_cast<int>(inequality_min_vector.size());

  std::cout << "minimize   1/2 x^T H x + x^T g" << std::endl;
  std::cout << "subject to lb <= x <= ub" << std::endl;
  std::cout << "           lbA <= Ax <= ubA" << std::endl;

  Eigen::SparseMatrix<double> hessian = eval_weight_matrix.sparseView();
  Eigen::VectorXd gradient = eval_coeff_vector;
  Eigen::MatrixXd linear_matrix_org(x_size + eq_size + ieq_size, x_size);
  linear_matrix_org << Eigen::MatrixXd::Identity(x_size, x_size), equality_matrix, inequality_matrix;
  Eigen::SparseMatrix<double> linear_matrix = linear_matrix_org.sparseView();
  Eigen::VectorXd lower_bound(x_size + eq_size + ieq_size);
  lower_bound << state_min_vector, equality_coeff_vector, inequality_min_vector;
  Eigen::VectorXd upper_bound(x_size + eq_size + ieq_size);
  upper_bound << state_max_vector, equality_coeff_vector, inequality_max_vector;

  OsqpEigen::Solver solver{};
  // solver.settings()->set<OsqpEigen::Setting>();
  solver.settings()->setWarmStart(true);
  solver.data()->setNumberOfVariables(x_size);
  solver.data()->setNumberOfConstraints(x_size + eq_size + ieq_size);
  if (not solver.data()->setHessianMatrix(hessian)) {
      std::cerr << "Error setting Hessian Matrix" << std::endl;
  }
  if (not solver.data()->setGradient(gradient)) {
      std::cerr << "Error setting Gradient" << std::endl;
  }
  if (not solver.data()->setLinearConstraintsMatrix(linear_matrix)) {
      std::cerr << "Error setting Linear Constraint Matrix" << std::endl;
  }
  if (not solver.data()->setLowerBound(lower_bound)) {
      std::cerr << "Error setting Lower Bound" << std::endl;
  }
  if (not solver.data()->setUpperBound(upper_bound)) {
      std::cerr << "Error setting Upper Bound" << std::endl;
  }
  if (not solver.initSolver()) {
      std::cerr << "Error init Solver" << std::endl;
  }

  if (not solver.solve()) {
      std::cerr << "Could not solve the QP problem." << std::endl;
      return std::numeric_limits<double>::infinity();
  }

  result_vector = solver.getSolution();

  return 0.0;
}
// }}}

int main(int argc, char** argv)
{
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();


  hrp::dvector state_min_vector(12);
  state_min_vector <<
-0.147462,-0.109559,-0.247033,-0.286287,-0.730731,-1.17518,-1.35934,-1.31123,-1.27514,-1.27784,-1.2743,-1.27316;

  hrp::dvector  state_max_vector(12);
  state_max_vector <<
1.41809,1.45599,1.31852,1.27926,0.83482,0.390375,0.206213,0.25432,0.290411,0.287713,0.291249,0.292391;

  hrp::dmatrix eval_weight_matrix(12, 12);
  eval_weight_matrix <<
1,0,0,0,0,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0,0,0,0,
0,0,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,1,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,0,
0,0,0,0,0,0,0,1,0,0,0,0,
0,0,0,0,0,0,0,0,1,0,0,0,
0,0,0,0,0,0,0,0,0,1,0,0,
0,0,0,0,0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,0,0,0,0,1;

  hrp::dvector eval_coeff_vector(12);
  eval_coeff_vector << 0,0,0,0,0,0,0,0,0,0,0,0;

  hrp::dmatrix equality_matrix(3, 12);
  equality_matrix <<
3.40342e-06,0.0215941,0.309683,0.528402,0.137753,0.00256421,0,0,0,0,0,0,
-0.000729304,-0.79947,-4.13328,1.81394,2.97347,0.146063,0,0,0,0,0,0,
0,0,0.000260417,0.0617188,0.438021,0.438021,0.0617187,0.000260417,0,0,0,0;

  hrp::dvector equality_coeff_vector(3);
  equality_coeff_vector <<
0,0,0.129184;

  hrp::dmatrix inequality_matrix(11, 12);
  inequality_matrix <<
-9,9,0,0,0,0,0,0,0,0,0,0,
0,-9,9,0,0,0,0,0,0,0,0,0,
0,0,-9,9,0,0,0,0,0,0,0,0,
0,0,0,-9,9,0,0,0,0,0,0,0,
0,0,0,0,-9,9,0,0,0,0,0,0,
0,0,0,0,0,-9,9,0,0,0,0,0,
0,0,0,0,0,0,-9,9,0,0,0,0,
0,0,0,0,0,0,0,-9,9,0,0,0,
0,0,0,0,0,0,0,0,-9,9,0,0,
0,0,0,0,0,0,0,0,0,-9,9,0,
0,0,0,0,0,0,0,0,0,0,-9,9;

  hrp::dvector inequality_min_vector(11);
  inequality_min_vector <<
-3.65888,-5.23726,-4.35329,-8,-8,-5.65746,-3.56704,-3.67518,-4.02428,-3.96818,-3.98973;

  hrp::dvector inequality_max_vector(11);
  inequality_max_vector <<
4.34112,2.76274,3.64671,3.28213e-06,-1.43026e-06,2.34254,4.43296,4.32482,3.97572,4.03182,4.01027;

  hrp::dvector tmp_dp_modified(12);

  double opt = solve_qp(state_min_vector, state_max_vector,
          eval_weight_matrix, eval_coeff_vector,
          equality_matrix, equality_coeff_vector,
          inequality_matrix, inequality_min_vector, inequality_max_vector,
          tmp_dp_modified);

  std::cout << "ans: " << tmp_dp_modified.transpose() << std::endl;
#warning always zero
  std::cout << "opt: " << opt << std::endl;

  std::cout << "time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count() * 1e-9 << std::endl;

  return 0;
}
