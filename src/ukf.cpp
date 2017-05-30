#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  //x_ = VectorXd(5);
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  //P_ = MatrixXd(5, 5);
  P_ = MatrixXd(VectorXd::Ones(n_x_).asDiagonal());

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.03;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO: Done

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  time_us_ = 0;
  is_initialized_ = false;

  // initialize process noise matrix
  Q_ = MatrixXd::Zero(2,2);
  Q_(0,0) = std_a_*std_a_;
  Q_(1,1) = std_yawdd_*std_yawdd_;

  // initialize radar measurement noise matrix
  R_= MatrixXd::Zero(3,3);
  R_(0,0) = std_radr_ * std_radr_;
  R_(1,1) = std_radphi_ * std_radphi_;
  R_(2,2) = std_radrd_ * std_radrd_;

  // initialize the predict sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_,2 * n_aug_+1);

  weights_ = VectorXd(2*n_aug_+1);
  double coeff = 0.5/(lambda_+n_aug_);
  weights_.fill(coeff);
  weights_[0] = lambda_/(lambda_+n_aug_);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO: Done

  Complete this function! Make sure you switch between lidar and radar
  measurements.
   */

  //record the timestamp and calculate the delta t
  double delta_t = (meas_package.timestamp_ - time_us_) / 1.0e6;
  time_us_ = meas_package.timestamp_;

  //for fist measurement, initialize the state x
  if (false == is_initialized_) {
    ProcessFirstMeasurement(meas_package);
    return;
  }else{

    //do prediction
    Prediction(delta_t);

    // do update
    switch (meas_package.sensor_type_) {
      case MeasurementPackage::LASER: {
        UpdateLidar(meas_package);
        break;
      }
      case MeasurementPackage::RADAR: {
        UpdateRadar(meas_package);
        break;
      }
      default: {
        cout << "Unexpected SensorType : " << meas_package.sensor_type_ << endl;
        return;
      }
    }
  }
}

/**
 * process first measurement
 * @param meas_package The latest measurement data of either radar or laser
 */
void UKF::ProcessFirstMeasurement(MeasurementPackage meas_package) {

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    /**
     Convert radar from polar to cartesian coordinates and initialize state.
     */
    if (meas_package.raw_measurements_.size() != 3) {
      cout << "Invalid RADA measurement data !!!" << endl;
      return;
    }

    double r = meas_package.raw_measurements_[0];
    double phi = meas_package.raw_measurements_[1];
    x_.fill(0);
    x_[0] = r * cos(phi);
    x_[1] = r * sin(phi);
    is_initialized_ = true;

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    /**
     Initialize state.
     */
    if (meas_package.raw_measurements_.size() != 2) {
      cout << "Invalid LIDA measurement data !!!" << endl;
      return;
    }
    //initialize x
    x_.fill(0);
    x_[0] = meas_package.raw_measurements_[0];
    x_[1] = meas_package.raw_measurements_[1];
    is_initialized_ = true;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:Done

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // generate sigma points
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  // calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  // calculate the square root of lambda + n_x_
  double coef = sqrt(3.0);  // sqrt(3)

  // assign each column with x_ value
  Xsig = x_ * MatrixXd::Ones(1, 2 * n_x_ + 1);
  Xsig.block(0, 1, n_x_, n_x_) += coef * A;
  Xsig.block(0, 1 + n_x_, n_x_, n_x_) -= coef * A;

  // create augmented mean vector
  VectorXd x_aug_mean = VectorXd::Zero(n_aug_);
  x_aug_mean.head(n_x_) = x_;
  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  // create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;

  // create square root matrix of P_aug
  MatrixXd L = P_aug.llt().matrixL();
  // assign each column with x_aug_mean
  Xsig_aug = x_aug_mean * MatrixXd::Ones(1, 2 * n_aug_ + 1);
  Xsig_aug.block(0,1,n_aug_,n_aug_) += coef * L;
  Xsig_aug.block(0,1 + n_aug_,n_aug_,n_aug_) -= coef * L;

  //predict sigma points

  double dt2= delta_t * delta_t;

  for(unsigned int i = 0; i< 2* n_aug_; i++){
    VectorXd col = Xsig_aug.col(i);
    double px = col[0];
    double py = col[1];
    double v = col[2];
    double yaw = col[3];
    double yawd = col[4];
    double nu_a = col[5];
    double nu_yawdd = col[6];
    if (fabs(yawd) > 0.001) {
      Xsig_pred_(0, i) = px + v * cos(yaw) * delta_t
          + 0.5 * dt2 * cos(yaw) * nu_a;
      Xsig_pred_(1, i) = py + v * sin(yaw) * delta_t
          + 0.5 * dt2 * sin(yaw) * nu_a;

    } else {
      Xsig_pred_(0, i) = px + v/ yawd*(sin(yaw + yawd*delta_t) - sin(yaw))
          + 0.5 * dt2 * cos(yaw) * nu_a;
      Xsig_pred_(1, i) = py + v/ yawd*(-cos(yaw + yawd*delta_t) + cos(yaw))
          + 0.5 * dt2 * sin(yaw) * nu_a;
    }

    Xsig_pred_(2,i) = v + delta_t*nu_a;
    Xsig_pred_(3,i) = yaw + yawd * nu_yawdd + 0.5 * dt2 * nu_yawdd;
    Xsig_pred_(4,i) = yawd + delta_t * nu_yawdd;
  }

  //predict state mean
  x_ = Xsig_pred_ * weights_;
  //predict state covariance matrix
  MatrixXd Diff = Xsig_pred_ - x_ * MatrixXd::Ones(1,2*n_aug_ +1);
  P_ = Diff * weights_.asDiagonal()* Diff.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;
  MatrixXd Zsig_pred = MatrixXd(n_z, 2*n_aug_ +1);
  //transform sigma points into measurement space
  for(unsigned int i = 0; i< 2* n_aug_ +1; i++){
    VectorXd col = Xsig_pred_.col(i);
    double px = col[0];
    double py = col[1];
    double v = col[2];
    double yaw = col[3];
    double yawd = col[4];

    double r, phi, r_dot;
    r = sqrt(px*px + py*py);
    phi = atan2(py,px);
    r_dot = (px * v*cos(yaw) + py*v*sin(yaw))/r;
    Zsig_pred.col(i)<< r,phi,r_dot;
  }
  VectorXd z_pred_mean = VectorXd::Zero(n_z);

  z_pred_mean = Zsig_pred * weights_;

  MatrixXd S = MatrixXd::Zero(n_z,n_z);




}
