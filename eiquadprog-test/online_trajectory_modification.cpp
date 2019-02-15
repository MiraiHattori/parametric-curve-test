#include "bspline.hpp"
#include "3rdparty/eiquadprog.hpp"

/*
hrp::dmatrix calcDpFromDr(const hrp::dvector& dr, const double& t_current, const double& t_hit, const std::vector<BSpline::BSpline>& bsplines, const std::vector<double>& p, const double& epsilon = 1.0e-6)
{
  int online_modified_min_id = -1;
  // bsplineのサイズは1以上，その関節は動き始める前提
  BSpline::BSpline bspline = bsplines.at(0);
  hrp::dvector coeff_vector = bspline.calcCoeffVector(t_current);
  int id_max = coeff_vector.size();
  for (int i = 0; i < coeff_vector.size(); i++)
  {
    if (std::abs(coeff_vector[i]) > epsilon)
    {
      online_modified_min_id = i;
      break;
    }
  }
  int online_modified_max_id_1 = -1;
  for (int i = 0; i < coeff_vector.size(); i++)
  {
    if (std::abs(coeff_vector[coeff_vector.size() - i]) > epsilon)
    {
      online_modified_min_id = coeff_vector.size() - i;
      break;
    }
  }
  // online_modified_max_id_1は実際のonline_modified_max_idより1大きい
  int c = online_modified_max_id_1 - online_modified_min_id;  // ここは+1する必要がない
  // online_modified_links ikを解く*limb*のlinkのリスト
  // online_modified_jlist online_modified_linksの:jointのリスト
  // k: online_modified_jlistの長さ
  // current_pose: t_hitでの関節角度+rootlink6自由度の計算, fix-leg-to-coordsなどをして現在姿勢を計算(ikの初期値の姿勢)
  // 元のEuslispコードではここで*robot*に関節角を設定している
  std::vector<double> current_pose{};
  for (size_t i = 0; i < bsplines.size(); i++)
  {
    BSpline::BSpline bspline = bsplines.at(i);
    hrp::dvector ps(id_max);
    for (int j = 0; j < id_max; j++)
    {
      ps[j] = p.at(id_max * i + j);
    }
#warning jointの数をちゃんと取得
    // joint数33, bsplinesのサイズが39(virtualjoint6dofがあるため)
    // TODO virtualjoint6dofをちゃんと設定する(ロボットを地面に接地させる)
    if (i < 33)
    {
      current_pose.push_back(bspline.calc(t_hit, ps));
    }
  }
  // dq
  // まずはラケット先端をendcoordsにする必要がある
  // ラケット先端をtargetに移動させる
  // その上でIKを解く
  // 関節角の差分をstd::vector<double>とかにして返す

  // dp 返り値の宣言
  std::vector<double> dp{p.size(), 0.0}; // id_max * ((length jlist) + 6) + 1(t_hit)
  // dp-modifiedを計算
  Eigen::VectorXd initial_state = Eigen::VectorXd::Zero(c);
  Eigen::MatrixXd initial_matrix = Eigen::MatrixXd::Zero(2, c);
}
*/

int main()
{
  return 0;
}
