#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
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
  VectorXd rsme(4);

  rsme << 0, 0, 0, 0;

  // Check availability of the the following inputs
  // Estimation vector size must not be zero

  float ti = estimations.size();

  if(ti == 0){
  	cout << "Error : Input Empty" << endl;
  	return rsme;
  }
  // The estimation vector size must equal ground truth vector size
  if(ti != ground_truth.size()){
  	cout << "Error : Size of estimation and ground_truth do not match" << endl;
  	return rsme;
  }

  // Squared Residual

  for(unsigned int i=0; i < ti; ++i){

  	VectorXd residual = estimations[i] - ground_truth[i];
  	residual = residual.array() * residual.array();
  	rsme = rsme + residual;
  }

  // Mean
  rsme = rsme/ti;

  // Squared root
  rsme = rsme.array().sqrt();

  return rsme;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

  MatrixXd Hj(3,4);


  // Recover State parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Precompute a set of terms to avoid repeated calculation
  float c1 = px * px + py * py;
  float c2 = sqrt(c1);
  float c3 = c1 * c2;

  // Check division by zero
  if(fabs(c1) < 0.0001){
  	cout << "CalculateJacobian() - Error - Div. by 0" << endl;
  	return Hj;
  }

  // Compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;

  
}
