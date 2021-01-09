#include "qpOASES.hpp"
#include <iostream>
#include <iomanip>
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

  std::cout << "minimize   1/2 x^T H x + x^T g" << std::endl;
  std::cout << "subject to lb <= x <= ub" << std::endl;
  std::cout << "           lbA <= Ax <= ubA" << std::endl;

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
  hrp::dmatrix A_Eigen(eq_size + ieq_size, x_size);
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
  hrp::dvector lbA_Eigen(eq_size + ieq_size);
  lbA_Eigen << equality_coeff_vector, inequality_min_vector;
  real_t* lbA = lbA_shptr.get();
  setUpQPParam(lbA, lbA_Eigen);

  boost::shared_ptr<real_t> ubA_shptr(new real_t[eq_size + ieq_size]);
  hrp::dvector ubA_Eigen(eq_size + ieq_size);
  ubA_Eigen << equality_coeff_vector, inequality_max_vector;
  real_t* ubA = ubA_shptr.get();
  setUpQPParam(ubA, ubA_Eigen);

  std::cout << "H" << std::endl;
  for (int i = 0; i < x_size; i++) {
      for (int j = 0; j < x_size; j++) {
          std::cout << H[i * x_size + j] << " ";
      }
      std::cout << std::endl;
  }
  std::cout << "g" << std::endl;
  for (int i = 0; i < x_size; i++) {
      std::cout << g[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "A" << std::endl;
  for (int i = 0; i < eq_size + ieq_size; i++) {
      for (int j = 0; j < x_size; j++) {
          std::cout << A[i * x_size + j] << " ";
      }
      std::cout << std::endl;
  }
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
  // options.printLevel = PL_NONE;
  // options.printLevel = PL_LOW;
  // options.printLevel = PL_MEDIUM;
  // options.printLevel = PL_HIGH;
  // options.printLevel = PL_TABULAR;
  options.printLevel = PL_DEBUG_ITER;
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

  std::cout << "lb:" << std::endl;
  for (int i = 0; i < state_min_vector.size(); i++) {
      std::cout << std::fixed << std::setprecision(6) << state_min_vector[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "x:" << std::endl;
  for (int i = 0; i < result_vector.size(); i++) {
      std::cout << std::fixed << std::setprecision(6) << result_vector[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "ub:" << std::endl;
  for (int i = 0; i < state_max_vector.size(); i++) {
      std::cout << std::fixed << std::setprecision(6) << state_max_vector[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "1/2 H^T x H + x^T g :" << std::endl;
  std::cout << std::fixed << std::setprecision(6) << result_vector.transpose() * eval_weight_matrix * result_vector / 2.0 + result_vector.transpose() * eval_coeff_vector << std::endl;
  std::cout << "lbA:" << std::endl;
  for (int i = 0; i < lbA_Eigen.size(); i++) {
      std::cout << std::fixed << std::setprecision(6) << lbA_Eigen[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Ax:" << std::endl;
  for (int i = 0; i < (A_Eigen * result_vector).size(); i++) {
      std::cout << std::fixed << std::setprecision(6) << (A_Eigen * result_vector)[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "ubA:" << std::endl;
  for (int i = 0; i < ubA_Eigen.size(); i++) {
      std::cout << std::fixed << std::setprecision(6) << ubA_Eigen[i] << " ";
  }
  std::cout << std::endl;

  return qp.getObjVal();
}
// }}}

int main(int argc, char** argv)
{
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  using namespace qpOASES;

  // https://hackmd.io/GjSw9VZ5TOyC5kfP5llMkQ?view

  /* Setup data of QP. */

  hrp::dvector state_min_vector(12);
  state_min_vector <<
-0.147462,-0.109559,-0.247033,-0.286287,-0.730731,-1.17518,-1.35934,-1.31123,-1.27514,-1.27784,-1.2743,-1.27316;

  hrp::dvector state_max_vector(12);
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

  std::cout << "opt: " << tmp_dp_modified.transpose() << std::endl;

  std::cout << "time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count() * 1e-9 << std::endl;

  return 0;
}
