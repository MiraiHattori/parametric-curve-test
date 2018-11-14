#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd mat(1, 10);
    mat << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::VectorXd vec(5);
    vec << 10, 11, 12, 13, 14;
    std::cout << mat.block(0, 0, 1, 5) << std::endl;
    mat.block(0, 0, 1, 5) = vec.transpose();
    std::cout << mat << std::endl;
    std::cout << mat.row(0).block(0, 0, 0, 9) << std::endl;
    Eigen::MatrixXd mat2(10, 1);
    mat2 << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
    std::cout << mat2.col(0) << std::endl;
    std::cout << mat2.col(0).segment(0, 10) << std::endl;
    Eigen::VectorXd vec2(10);
    vec2 = mat2.col(0);
    std::cout << vec2 << std::endl;
    std::cout << vec2.segment(0, 10) << std::endl;
    return 0;
}
