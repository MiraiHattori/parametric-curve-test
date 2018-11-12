#include "bspline.hpp"
#include <iostream>
#include <random>

int main()
{
    // std::cout << BSpline::calcBSplineCoeffVector(0.0, 3, 3) << std::endl;

    int id_max = 8;
    // recursive_order, recursive_cnt, id, id_max, t_min, t_max
    BSpline::BSpline a(3, 3, id_max, id_max, 0.0, 1.0);
    hrp::dvector gain = hrp::dvector::Zero(id_max);
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_real_distribution<> rand(0.0, 1.0);
    for (int i = 0; i < id_max; i++) {
        gain[i] = 2.0 * rand(mt) - 1.0;
    }
    a.varifyCalcDelta(gain);
    std::cout << a.calcDeltaMatrixForKeepRecursiveOrder(1) << std::endl;
    return 0;
}
