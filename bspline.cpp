#include "bspline.hpp"
#include <iostream>
#include <random>

int main()
{
    int id_max = 8;
    double t_min = 0.0;
    double t_max = 1.0;
    // recursive_order, recursive_cnt, id, id_max, t_min, t_max
    BSpline::BSpline a(3, 3, id_max, id_max, t_min, t_max);
    hrp::dvector gain = hrp::dvector::Zero(id_max);
    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_real_distribution<> rand(0.0, 1.0);
    for (int i = 0; i < id_max; i++) {
        gain[i] = 2.0 * rand(mt) - 1.0;
    }
    a.varifyCalcDelta(gain);
    std::cout << a.calcDeltaMatrixForKeepRecursiveOrder(1) << std::endl;

    double tm = rand(mt);
    int n = 1;
    hrp::dvector ref_dcv = a.calcDeltaCoeffVector(tm, n);
    hrp::dmatrix m_j = a.calcCoeffMatrixForTimeVector(tm);
    hrp::dmatrix m = BSpline::calcDeltaCoeffMatrixForTimeVector(a.recursiveOrder());
    hrp::dmatrix ret = m;
    for (int i = 0; i < n - 1; i++) {
        ret = m * ret;
    }
    hrp::dmatrix d_n = ret;
    hrp::dvector t = BSpline::calcTVector(tm, true, hrp::dvector::Zero(4), a.recursiveOrder());
    hrp::dvector dcv = (m_j * d_n) * t;
    std::cout << 1 << std::endl;
    std::cout << m_j << std::endl;
    std::cout << d_n << std::endl;
    std::cout << d_n * t << std::endl;
    std::cout << ref_dcv.transpose() << " - " <<  dcv.transpose() << " = " << (ref_dcv - dcv).transpose() << "(" << (ref_dcv - dcv).norm() << ")" << std::endl;
    std::cout << (ref_dcv - dcv).norm() << std::endl;
    return 0;
}
