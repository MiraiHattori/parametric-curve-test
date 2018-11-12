#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(0, 3);
    Eigen::VectorXd a0(3);
    Eigen::VectorXd a1(3);
    Eigen::VectorXd a2(3);
    a0 << 1, 2, 3;
    a1 << 4, 5, 6;
    a2 << 7, 8, 9;
    for (const Eigen::VectorXd& elem : {a0, a1, a2}) {
        // std::cout << elem.transpose() << std::endl;
        Eigen::MatrixXd b = Eigen::MatrixXd::Zero(a.rows() + 1, a.cols());
        for (int i = 0; i < a.rows() + 1; i++) {
            if (i == a.rows()) {
                b.row(i) = elem;
            } else {
                b.row(i) = a.row(i);
            }
        }
        a = b;
    }
    std::cout << a << std::endl;
}
