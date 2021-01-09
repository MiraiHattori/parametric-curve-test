#include "qpOASES.hpp"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
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

// {{{ qpOASESWrapper
void setUpQPParam(qpOASES::real_t* raw_mat, const hrp::dmatrix& hrp_dmat)
{
    for (int i = 0; i < hrp_dmat.rows(); i++) {
        for (int j = 0; j < hrp_dmat.cols(); j++) {
            raw_mat[i * hrp_dmat.cols() + j] = hrp_dmat(i, j);
        }
    }
}

void setUpQPParam(qpOASES::real_t* raw_arr, const hrp::dvector& hrp_dvec)
{
    for (int i = 0; i < hrp_dvec.size(); i++) {
        raw_arr[i] = hrp_dvec[i];
    }
}

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
  int x_size = eval_weight_matrix.rows();
  int eq_size = equality_coeff_vector.size();
  int ieq_size = inequality_min_vector.size();
  std::cout << "x_size: " << x_size << std::endl;
  std::cout << "eq_size: " << eq_size << std::endl;
  std::cout << "ieq_size: " << ieq_size << std::endl;
  using namespace qpOASES;
  boost::shared_ptr<real_t> H_shptr(new real_t[x_size * x_size]);
  real_t* H = H_shptr.get();
  setUpQPParam(H, eval_weight_matrix);

  boost::shared_ptr<real_t> g_shptr(new real_t[x_size]);
  real_t* g = g_shptr.get();
  setUpQPParam(g, eval_coeff_vector);

  boost::shared_ptr<real_t> A_shptr(new real_t[(eq_size + ieq_size) * x_size]);
  Eigen::MatrixXd A_Eigen(eq_size + ieq_size, x_size);
  A_Eigen << equality_matrix, inequality_matrix;
  real_t* A = A_shptr.get();
  setUpQPParam(A, A_Eigen);

  boost::shared_ptr<real_t> lb_shptr(new real_t[x_size]);
  real_t* lb = lb_shptr.get();
  setUpQPParam(lb, state_min_vector);

  boost::shared_ptr<real_t> ub_shptr(new real_t[x_size]);
  real_t* ub = ub_shptr.get();
  setUpQPParam(ub, state_max_vector);

  boost::shared_ptr<real_t> lbA_shptr(new real_t[eq_size + ieq_size]);
  Eigen::VectorXd lbA_Eigen(eq_size + ieq_size);
  lbA_Eigen << equality_coeff_vector, inequality_min_vector;
  real_t* lbA = lbA_shptr.get();
  setUpQPParam(lbA, lbA_Eigen);

  boost::shared_ptr<real_t> ubA_shptr(new real_t[eq_size + ieq_size]);
  Eigen::VectorXd ubA_Eigen(eq_size + ieq_size);
  ubA_Eigen << equality_coeff_vector, inequality_max_vector;
  real_t* ubA = ubA_shptr.get();
  setUpQPParam(ubA, ubA_Eigen);

  std::cout << "H" << std::endl;
  for (int i = 0; i < x_size * x_size; i++) {
      std::cout << H[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "g" << std::endl;
  for (int i = 0; i < x_size; i++) {
      std::cout << g[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "A" << std::endl;
  for (int i = 0; i < (eq_size + ieq_size) * x_size; i++) {
      std::cout << A[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "lb" << std::endl;
  for (int i = 0; i < x_size; i++) {
      std::cout << lb[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "ub" << std::endl;
  for (int i = 0; i < x_size; i++) {
      std::cout << ub[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "lbA" << std::endl;
  for (int i = 0; i < eq_size + ieq_size; i++) {
      std::cout << lbA[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "ubA" << std::endl;
  for (int i = 0; i < eq_size + ieq_size; i++) {
      std::cout << ubA[i] << " ";
  }
  std::cout << std::endl;

  QProblem qp(x_size, eq_size + ieq_size);
  Options options;
  options.printLevel = PL_NONE;
  qp.setOptions(options);
  sparse_int_t nWSR = 100;
  real_t max_cputime_sec[1] = { 1e-3 };
  returnValue qp_init = qp.init(H, g, A, lb, ub, lbA, ubA, nWSR, max_cputime_sec);
  switch (qp_init) {
  case SUCCESSFUL_RETURN:
      std::cout << "SUCCESSFUL_RETURN" << std::endl;
      break;
  case RET_INIT_SUCCESSFUL:
      std::cout << "RET_INIT_SUCCESSFUL" << std::endl;
      break;
  case RET_INIT_FAILED:
      std::cout << "RET_INIT_FAILED" << std::endl;
      break;
  case RET_INIT_FAILED_TQ:
      std::cout << "RET_INIT_FAILED_TQ" << std::endl;
      break;
  case RET_INIT_FAILED_CHOLESKY:
      std::cout << "RET_INIT_FAILED_CHOLESKY" << std::endl;
      break;
  case RET_INIT_FAILED_HOTSTART:
      std::cout << "RET_INIT_FAILED_HOTSTART" << std::endl;
      break;
  case RET_INIT_FAILED_INFEASIBILITY:
      std::cout << "RET_INIT_FAILED_INFEASIBILITY" << std::endl;
      break;
  case RET_INIT_FAILED_UNBOUNDEDNESS:
      std::cout << "RET_INIT_FAILED_UNBOUNDEDNESS" << std::endl;
      break;
  case RET_INIT_FAILED_REGULARISATION:
      std::cout << "RET_INIT_FAILED_REGULARISATION" << std::endl;
      break;
  default:
      std::cout << "unknown qp init status: " << qp_init << std::endl;
      break;
  };
  boost::shared_ptr<real_t> x_opt_shptr(new real_t[x_size]);
  real_t* x_opt = x_opt_shptr.get();
  returnValue primal_solution_value = qp.getPrimalSolution(x_opt);
  if (primal_solution_value == SUCCESSFUL_RETURN) {
      std::cout << "SUCCESSFUL_RETURN" << std::endl;
  } else if (primal_solution_value == RET_QP_NOT_SOLVED) {
      std::cout << "RET_QP_NOT_SOLVED" << std::endl;
  } else {
    std::cout << "unknown qp primal solution status: " << primal_solution_value << std::endl;
  }
  for (int i = 0; i < result_vector.size(); i++) {
      result_vector[i] = x_opt[i];
  }
  return qp.getObjVal();
}
// }}}

int main(int argc, char** argv)
{
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  using namespace qpOASES;

  // https://hackmd.io/GjSw9VZ5TOyC5kfP5llMkQ?view

  /* Setup data of QP. */

  hrp::dvector state_min_vector(2);
  state_min_vector <<
0.5, -2.0;

  hrp::dvector  state_max_vector(2);
  state_max_vector <<
5.0, 2.0;

  hrp::dmatrix eval_weight_matrix(2, 2);
  eval_weight_matrix <<
1.0, 0.0,
0.0, 0.5;

  hrp::dvector eval_coeff_vector(2);
  eval_coeff_vector << 1.5, 1.0;

  hrp::dmatrix equality_matrix = hrp::dmatrix::Zero(1, 2);

  hrp::dvector equality_coeff_vector = hrp::dvector::Zero(1);

  hrp::dmatrix inequality_matrix(1, 2);
  inequality_matrix <<
1.0, 1.0;

  hrp::dvector inequality_min_vector(1);
  inequality_min_vector <<
-1.0;

  hrp::dvector inequality_max_vector(1);
  inequality_max_vector <<
2.0;

  hrp::dvector tmp_dp_modified(2);

  double opt = solve_qp(state_min_vector, state_max_vector,
          eval_weight_matrix, eval_coeff_vector,
          equality_matrix, equality_coeff_vector,
          inequality_matrix, inequality_min_vector, inequality_max_vector,
          tmp_dp_modified);

  std::cout << "opt: " << tmp_dp_modified.transpose() << std::endl;

  std::cout << "time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count() * 1e-9 << std::endl;

  return 0;
}
