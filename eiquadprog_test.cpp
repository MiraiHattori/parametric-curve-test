#include "3rdparty/eiquadprog.hpp"
#include <iostream>

/*
 * minimize 1/2 x^T G x + g0^T x
 * subject to CE^T x + ce0 = 0
 * subject to CI^T x + ci0 >= 0
 */
int main()
{

  /* Setup data of QP. */
  Eigen::MatrixXd G(2, 2);
  G << 4.0, 1.0, 1.0, 2.0; // 正定値ヘッセ行列
  Eigen::VectorXd g0(2); // 勾配ベクトル(?)
  g0 << 1.0, 1.0;
  Eigen::MatrixXd CE(2, 1); // constraint matrix
  CE << 1.0, 1.0;
  Eigen::VectorXd ce0(1); // 等式
  ce0 << -1.0;
  // dummies
  Eigen::MatrixXd CI(2, 1); // constraint matrix
  CI << 1.0, 1.0;
  Eigen::VectorXd ci0(1); // 等式
  ci0 << 1.0;

  Eigen::VectorXd x(2);
  std::cout << Eigen::solve_quadprog(G, g0, CE, ce0, CI, ci0, x) << std::endl;
  std::cout << x.transpose() << std::endl;

  return 0;
}
