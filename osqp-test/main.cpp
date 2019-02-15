#include <iostream>
#include <OsqpEigen/OsqpEigen.h>

int main()
{

    // allocate QP problem matrices and vectores
    Eigen::MatrixXd hessian_org(12, 12);
    hessian_org <<
1,0,0,0,0,0,0,0,0,0,0,0,
0,1.1,0,0,0,0,0,0,0,0,0,0,
0,0,1.2,0,0,0,0,0,0,0,0,0,
0,0,0,1.3,0,0,0,0,0,0,0,0,
0,0,0,0,1.4,0,0,0,0,0,0,0,
0,0,0,0,0,1.5,0,0,0,0,0,0,
0,0,0,0,0,0,1.6,0,0,0,0,0,
0,0,0,0,0,0,0,1.7,0,0,0,0,
0,0,0,0,0,0,0,0,1.8,0,0,0,
0,0,0,0,0,0,0,0,0,1.9,0,0,
0,0,0,0,0,0,0,0,0,0,2,0,
0,0,0,0,0,0,0,0,0,0,0,2.1;
    Eigen::SparseMatrix<double> hessian = hessian_org.sparseView();

    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(12);

    Eigen::MatrixXd linear_matrix_org(26, 12);
    linear_matrix_org <<
3.40342e-06,0.0215941,0.309683,0.528402,0.137753,0.00256421,0,0,0,0,0,0,
-0.000729304,-0.79947,-4.13328,1.81394,2.97347,0.146063,0,0,0,0,0,0,
0,0,0.000260417,0.0617188,0.438021,0.438021,0.0617187,0.000260417,0,0,0,0,
1,0,0,0,0,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0,0,0,0,
0,0,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,1,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,0,
0,0,0,0,0,0,0,1,0,0,0,0,
0,0,0,0,0,0,0,0,1,0,0,0,
0,0,0,0,0,0,0,0,0,1,0,0,
0,0,0,0,0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,0,0,0,0,1,
-9,9,0,0,0,0,0,0,0,0,0,0,
0,-9,9,0,0,0,0,0,0,0,0,0,
0,0,-9,9,0,0,0,0,0,0,0,0,
0,0,0,-9,9,0,0,0,0,0,0,0,
0,0,0,0,-9,9,0,0,0,0,0,0,
0,0,0,0,0,-9,9,0,0,0,0,0,
0,0,0,0,0,0,-9,9,0,0,0,0,
0,0,0,0,0,0,0,-9,9,0,0,0,
0,0,0,0,0,0,0,0,-9,9,0,0,
0,0,0,0,0,0,0,0,0,-9,9,0,
0,0,0,0,0,0,0,0,0,0,-9,9;
    Eigen::SparseMatrix<double> linear_matrix = linear_matrix_org.sparseView();

    Eigen::VectorXd lower_bound(26);
lower_bound <<
-0,-0,-0.129184,
0.147462,0.109559,0.247033,0.286287,0.730731,1.17518,1.35934,1.31123,1.27514,1.27784,1.2743,1.27316,
3.65888,5.23726,4.35329,8,8,5.65746,3.56704,3.67518,4.02428,3.96818,3.98973;
    Eigen::VectorXd upper_bound(26);
upper_bound <<
-0,-0,-0.129184,
1.41809,1.45599,1.31852,1.27926,0.83482,0.390375,0.206213,0.25432,0.290411,0.287713,0.291249,0.292391,
4.34112,2.76274,3.64671,3.28213e-06,-1.43026e-06,2.34254,4.43296,4.32482,3.97572,4.03182,4.01027;


    OsqpEigen::Solver solver{};
    // solver.settings()->set<OsqpEigen::Setting>();
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(12);
    solver.data()->setNumberOfConstraints(26);
    if (not solver.data()->setHessianMatrix(hessian)) {
        std::cerr << "Error setting Hessian Matrix" << std::endl;
    }
    if (not solver.data()->setGradient(gradient)) {
        std::cerr << "Error setting Gradient" << std::endl;
    }
    if (not solver.data()->setLinearConstraintsMatrix(linear_matrix)) {
        std::cerr << "Error setting Linear Constraint Matrix" << std::endl;
    }
    if (not solver.data()->setLowerBound(lower_bound)) {
        std::cerr << "Error setting Lower Bound" << std::endl;
    }
    if (not solver.data()->setUpperBound(upper_bound)) {
        std::cerr << "Error setting Upper Bound" << std::endl;
    }
    if (not solver.initSolver()) {
        std::cerr << "Error setting Upper Bound" << std::endl;
    }
    std::cerr << "po" << std::endl;

    if (not solver.solve()) {
        std::cerr << "Could not solve the QP problem." << std::endl;
    }

    Eigen::VectorXd result = solver.getSolution();
    std::cout << result.transpose() << std::endl;

    return 0;
}
