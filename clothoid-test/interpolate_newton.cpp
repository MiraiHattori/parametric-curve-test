/*********************************************
 * ニュートン補間
 *********************************************/
#include <iostream>
#include <vector>

/*
 * 計算クラス
 */
class Calc
{
public:
    explicit Calc(const std::vector<double>& xs, const std::vector<double>& ys);
    void print();
    double interpolateNewton(const double& t);
    void addNewPoint(const double& x, const double& y);

private:
    std::vector<double> ws;
    std::vector<double> cs;
};

Calc::calc(const std::vector<double>& xs, const std::vector<double>& ys) {
    cs = std::vector<double>(xs.size());
    ws = std::vector<double>(xs.size());
}

void Calc::print()
{
    for (double t = xs[0]; t <= xs[xs.size() - 1]; t += 0.5) {
        std::cout << t << " " << interpolateNewton(t) << std::endl;
    }
}

double Calc::interpolateNewton(const double& t)
{
    std::vector<double> c(xs.size());
    std::vector<double> w(xs.size());
    // 係数
    for (int i = 0; i < xs.size(); i++) {
        w[i] = ys[i];
        for (int j = i - 1; j >= 0; j--) {
            w[j] = (w[j+1] - w[j]) / (xs[i] - xs[j]);
        }
        c[i] = w[0];
    }

    // 総和
    double sum = c[xs.size()-1];
    for (int i = xs.size() - 2; i >= 0; i--) {
        sum = sum * (t - xs[i]) + c[i];
    }

    return sum;
}

void addNewPoint(const double& xs, const double& ys) {
}

int main()
{
    try {
        std::vector<double> xs({0.0, 2.0, 3.0, 5.0, 8.0});
        std::vector<double> ys({0.8, 3.2, 2.8, 4.5, 1.9});
        Calc objCalc(xs, ys);

        objCalc.print();
    } catch (...) {
        std::cout << "例外発生！" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
