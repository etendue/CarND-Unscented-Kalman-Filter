#include <iostream>
#include "tools.h"
#include <assert.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   TODO: Done
   * Calculate the RMSE here.
   */

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size


  assert(estimations.size() != 0);
  assert (estimations.size() == ground_truth.size());

  unsigned int size = estimations.size();
  unsigned int dim = estimations[0].size();
  VectorXd rmse(dim);
  rmse.fill(0);

  for (unsigned int i = 0; i < size; i++) {
    assert(estimations[i].size()==ground_truth[i].size());
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  //calculate the mean
  rmse = rmse / size;

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}
