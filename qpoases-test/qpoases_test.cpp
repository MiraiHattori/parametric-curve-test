#include "qpOASES.hpp"
#include <iostream>

/*
 * minimize 1/2 x^T H x + c^T x
 * subject to Ax <= b
 */
int main()
{
  using namespace qpOASES;

  // https://hackmd.io/GjSw9VZ5TOyC5kfP5llMkQ?view

  /* Setup data of QP. */
  real_t H[2 * 2] = { 4.0, 1.0, 1.0, 2.0 };  // 正定値ヘッセ行列
  real_t A[1 * 2] = { 1.0, 1.0 };            // constraint matrix
  real_t g[2] = { 1.0, 1.0 };                // 勾配ベクトル
  real_t lb[2] = { 0.0, 0.0 };               // x lower bound
  real_t ub[2] = { 10000.0, 10000.0 };       // x upper bound
  real_t lbAx[1] = { 1.0 };                  // Ax lower bound; これとubAを等しくすることで等式条件となる
  real_t ubAx[1] = { 1.0 };

  /* Setting up QProblem object. */
  QProblem example(2, 1);

  Options options;
  example.setOptions(options);

  /* Solve first QP. */
  sparse_int_t nWSR = 10;
  example.init(H, g, A, lb, ub, lbAx, ubAx, nWSR);

  /* Get and print solution of first QP. */
  real_t xOpt[2];
  real_t yOpt[2 + 1];
  example.getPrimalSolution(xOpt);
  example.getDualSolution(yOpt);
  std::cout << std::endl;
  std::cout << "xOpt = [ " << xOpt[0] << ", " << xOpt[1] << " ];  yOpt = [ " << yOpt[0] << ", " << yOpt[1] << ", " << yOpt[2] << " ];  objVal = " << example.getObjVal() << std::endl;
  std::cout << std::endl;

  example.printOptions();
  example.printProperties();

  return 0;
}
