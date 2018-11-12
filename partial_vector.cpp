#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd mat(1, 10);
    mat << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::VectorXd vec(5);
    vec << 10, 11, 12, 13, 14;
    mat.block(0, 0, 1, 5) = vec.transpose();
    std::cout << mat << std::endl;
    return 0;
}
