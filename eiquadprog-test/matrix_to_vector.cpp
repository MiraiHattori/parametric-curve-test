#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd ret(4, 4);
    ret << 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20;
    Eigen::MatrixXd a(2, 2);
    a << 0, 1, 2, 3;
    a.transposeInPlace();
    ret.row(0) = Eigen::Map<Eigen::VectorXd>(a.data(), a.cols() * a.rows());
    std::cout << Eigen::Map<Eigen::VectorXd>(a.data(), a.cols() * a.rows()) << std::endl;
    std::cout << ret << std::endl;
    return 0;
}
