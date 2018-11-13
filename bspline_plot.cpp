#include "bspline.hpp"

int main()
{
  int id_max = 10;
  int recursive_order = 4;
  int recursive_cnt = 4;
  double t_min = 0.0;
  double t_max = 1.0;
  int plot_cnt = 100;
  BSpline::BSpline bspline(recursive_order, recursive_cnt, id_max, id_max, t_min, t_max);
  for (int i = 0; i < plot_cnt; i++)
  {
    double t = t_min + i * 1.0 / plot_cnt * (t_max - t_min);
    for (int j = 0; j < bspline.calcCoeffVector(t).size(); j++)
    {
      std::cout << j << " " << t << " " << bspline.calcCoeffVector(t)[j] << std::endl;
    }
  }
  return 0;
}
