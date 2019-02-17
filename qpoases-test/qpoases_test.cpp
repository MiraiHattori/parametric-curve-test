#include "qpOASES.hpp"
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <chrono>

/*
 * minimize 1/2 x^T H x + x^T g
 * subject to lb <= x <= ub
 *            lba <= Ax <= ubA
 */

void setUpQPParam(qpOASES::real_t* raw_mat, const Eigen::MatrixXd& eigen_mat)
{
    for (int i = 0; i < eigen_mat.rows(); i++) {
        for (int j = 0; j < eigen_mat.cols(); j++) {
            raw_mat[i * eigen_mat.cols() + j] = eigen_mat(i, j);
        }
    }
}

void setUpQPParam(qpOASES::real_t* raw_arr, const Eigen::VectorXd& eigen_arr)
{
    for (int i = 0; i < eigen_arr.size(); i++) {
        raw_arr[i] = eigen_arr[i];
    }
}

int main(int argc, char** argv)
{
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  using namespace qpOASES;

  // https://hackmd.io/GjSw9VZ5TOyC5kfP5llMkQ?view

  /* Setup data of QP. */
  boost::shared_ptr<real_t> H_shptr(new real_t[12 * 12]);
  Eigen::MatrixXd H_Eigen(12, 12);
  H_Eigen <<
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
  real_t* H = H_shptr.get();
  setUpQPParam(H, H_Eigen); // 正定値ヘッセ行列

  boost::shared_ptr<real_t> g_shptr(new real_t[12]);
  Eigen::VectorXd g_Eigen(12);
  g_Eigen << 0,0,0,0,0,0,0,0,0,0,0,0;
  real_t* g = g_shptr.get();
  setUpQPParam(g, g_Eigen); // 勾配ベクトル

  boost::shared_ptr<real_t> A_shptr(new real_t[14 * 12]);
  Eigen::MatrixXd A_Eigen(14, 12);
  A_Eigen <<
3.40342e-06,0.0215941,0.309683,0.528402,0.137753,0.00256421,0,0,0,0,0,0,
-0.000729304,-0.79947,-4.13328,1.81394,2.97347,0.146063,0,0,0,0,0,0,
0,0,0.000260417,0.0617188,0.438021,0.438021,0.0617187,0.000260417,0,0,0,0,
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
  real_t* A = A_shptr.get();
  setUpQPParam(A, A_Eigen);

  boost::shared_ptr<real_t> lb_shptr(new real_t[12]);
  Eigen::VectorXd lb_Eigen(12);
  lb_Eigen <<
-0.147462,-0.109559,-0.247033,-0.286287,-0.730731,-1.17518,-1.35934,-1.31123,-1.27514,-1.27784,-1.2743,-1.27316;
  real_t* lb = lb_shptr.get();
  setUpQPParam(lb, lb_Eigen);

  boost::shared_ptr<real_t> ub_shptr(new real_t[12]);
  Eigen::VectorXd ub_Eigen(12);
  ub_Eigen <<
1.41809,1.45599,1.31852,1.27926,0.83482,0.390375,0.206213,0.25432,0.290411,0.287713,0.291249,0.292391;
  real_t* ub = ub_shptr.get();
  setUpQPParam(ub, ub_Eigen);

  boost::shared_ptr<real_t> lbAx_shptr(new real_t[14]);
  Eigen::VectorXd lbAx_Eigen(14);
  lbAx_Eigen <<
0,0,0.129184,
-3.65888,-5.23726,-4.35329,-8,-8,-5.65746,-3.56704,-3.67518,-4.02428,-3.96818,-3.98973;
  real_t* lbAx = lbAx_shptr.get();
  setUpQPParam(lbAx, lbAx_Eigen);

  boost::shared_ptr<real_t> ubAx_shptr(new real_t[14]);
  Eigen::VectorXd ubAx_Eigen(14);
  ubAx_Eigen <<
0,0,0.129184,
4.34112,2.76274,3.64671,3.28213e-06,-1.43026e-06,2.34254,4.43296,4.32482,3.97572,4.03182,4.01027;
  real_t* ubAx = ubAx_shptr.get();
  setUpQPParam(ubAx, ubAx_Eigen);

  /* Setting up QProblem object. */
  QProblem example(12, 14);

  Options options;
  options.printLevel = PL_NONE;
  example.setOptions(options);

  /* Solve first QP. */
  sparse_int_t nWSR = 100;
  real_t max_cputime_sec[1] = { 1e-3 };
  returnValue qp_return_value = example.init(H, g, A, lb, ub, lbAx, ubAx, nWSR, max_cputime_sec);
  if (qp_return_value == SUCCESSFUL_RETURN) {
    std::cout << "qp succeeded" << std::endl;
  } else if (qp_return_value == RET_INIT_FAILED) {
    std::cout << "qp init failed" << std::endl;
  } else if (qp_return_value == RET_INIT_FAILED_CHOLESKY) {
    std::cout << "qp init failed cholesky" << std::endl;
  } else if (qp_return_value == RET_INIT_FAILED_TQ) {
    std::cout << "qp init failed tq" << std::endl;
  } else if (qp_return_value == RET_INIT_FAILED_HOTSTART) {
    std::cout << "qp init failed hotstart" << std::endl;
  } else if (qp_return_value == RET_INIT_FAILED_UNBOUNDEDNESS) {
    std::cout << "qp init failed unboundness" << std::endl;
  } else if (qp_return_value == RET_MAX_NWSR_REACHED){
    std::cout << "qp max nwsr reached" << std::endl;
  } else if (qp_return_value == RET_INVALID_ARGUMENTS) {
    std::cout << "qp init invalid arguments" << std::endl;
  } else {
    std::cout << "unknown qp status: " << qp_return_value << std::endl;
  }

  /* Get and print solution of first QP. */
  boost::shared_ptr<real_t> xOpt_shptr(new real_t[12]);
  real_t* xOpt = xOpt_shptr.get();
  returnValue primal_solution_value = example.getPrimalSolution(xOpt);
  if (primal_solution_value == SUCCESSFUL_RETURN) {
      std::cout << "qp primal solution succeeded" << std::endl;
  } else if (primal_solution_value == RET_QP_NOT_SOLVED) {
      std::cout << "qp not solved" << std::endl;
  } else {
    std::cout << "unknown qp primal solution status: " << primal_solution_value << std::endl;
  }

  boost::shared_ptr<real_t> yOpt_shptr(new real_t[12 + 14]);
  real_t* yOpt = yOpt_shptr.get();
  returnValue dual_solution_value = example.getDualSolution(yOpt);
  if (dual_solution_value == SUCCESSFUL_RETURN) {
      std::cout << "qp dual solution succeeded" << std::endl;
  } else if (dual_solution_value == RET_QP_NOT_SOLVED) {
      std::cout << "qp not solved" << std::endl;
  } else {
    std::cout << "unknown qp dual solution status: " << dual_solution_value << std::endl;
  }

  for (int i = 0; i < 12; i++) {
      std::cout << xOpt[i] << " ";
  }
  std::cout << std::endl;
  for (int i = 0; i < 12 + 14; i++) {
      std::cout << yOpt[i] << " ";
  }
  std::cout << std::endl;
  std::cout << example.getObjVal() << std::endl;

  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count() * 1e-9 << std::endl;

  //example.printOptions();
  //example.printProperties();

  return 0;
}
