#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    // TODO: Calculate the RMSE here.

    VectorXd rmse(ground_truth[0].size());
    rmse << VectorXd::Zero(ground_truth[0].size());

    //  the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size()) {
        cout << "Error. The estimation vector must be the same as the ground truth vector size." << endl;
        return rmse;
    }

    //  the estimation vector size should not be zero
    if (estimations.size() == 0) {
        cout << "Error. The estimation vector cannot have zero size." << endl;
        return rmse;
    }

    VectorXd residuals(estimations[0].size());

    //accumulate squared residuals

    for(int i=0; i < estimations.size(); ++i){
        residuals = estimations[i] - ground_truth[i];
        residuals = residuals.array() * residuals.array();
        rmse += residuals;

    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    // TODO: Calculate a Jacobian here.

    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //check division by zero
    if (px == 0 && py == 0) {
        cout << "Error. Division by zero." << endl;
        return Hj;
    }
    else {
        //compute the Jacobian matrix
        float px_py_sum = px*px + py*py;
        float px_py_sqrt = sqrt(px_py_sum);
        float px_py_cub_sqrt = px_py_sqrt * px_py_sqrt * px_py_sqrt;
        Hj << px / px_py_sqrt, py / px_py_sqrt, 0, 0,
                -py / px_py_sum, px / px_py_sum, 0, 0,
                (py * (vx*py - vy*px) / px_py_cub_sqrt), (px * (vy*px - vx*py) / px_py_cub_sqrt), (px / px_py_sqrt), (py / px_py_sqrt);

    }

    return Hj;

}
