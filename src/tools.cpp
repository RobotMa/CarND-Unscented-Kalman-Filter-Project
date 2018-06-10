#include <iostream>
#include "tools.h"

using Eigen::Vector4d;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    Vector4d vecRMSE(0, 0, 0, 0);
    if (estimations.size() != ground_truth.size() || estimations.size() < 1 )
    {
        return vecRMSE;
    }

    int n1 = estimations.size();

    Vector4d residual;

    for(int i = 0; i < n1; ++i)
    {
        residual = estimations.at(i) - ground_truth.at(i);
        residual = residual.array()*residual.array();
        vecRMSE += residual;
    }

    vecRMSE /= n1;
    vecRMSE = vecRMSE.array().sqrt();

    std::cout << "Final residual is " << vecRMSE << std::endl;

    return vecRMSE;
}
