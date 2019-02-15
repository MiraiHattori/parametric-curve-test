#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd m1(2, 2);
    m1 << 0, 1, 2, 3;
    Eigen::MatrixXd m2(2, 2);
    m2 << 4, 5, 6, 7;
    std::cout << m1 * m2 << std::endl;
    std::cout << m2 * m1 << std::endl;
    m1 *= m2;
    std::cout << m1 << std::endl;
    return 0;
}
