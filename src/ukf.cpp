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

  // initialize these value firsts
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  time_us_ = 0;
  is_initialized_ = false;
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
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.1;

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


  // initialize process noise matrix
  Q_ = MatrixXd::Zero(2,2);
  Q_(0,0) = std_a_*std_a_;
  Q_(1,1) = std_yawdd_*std_yawdd_;


  // initialize radar measurement noise matrix
  R_lidar_ = MatrixXd::Zero(2, 2);
  R_lidar_(0, 0) = std_laspx_ * std_laspx_;
  R_lidar_(1, 1) = std_laspy_ * std_laspy_;


  // initialize radar measurement noise matrix
  R_radar_= MatrixXd::Zero(3,3);
  R_radar_(0,0) = std_radr_ * std_radr_;
  R_radar_(1,1) = std_radphi_ * std_radphi_;
  R_radar_(2,2) = std_radrd_ * std_radrd_;

  // initialize the predict sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_,2 * n_aug_+1);

  weights_ = VectorXd(2*n_aug_+1);
  double coeff = 0.5/(lambda_+n_aug_);
  weights_.fill(coeff);
  weights_[0] = lambda_/(lambda_+n_aug_);

  ofs_nis_radar_.open("nis_radar.txt",ofstream::out);
  ofs_nis_lidar_.open("nis_lidar.txt",ofstream::out);

}

UKF::~UKF() {
  ofs_nis_radar_.close();
  ofs_nis_lidar_.close();
}

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
        if(use_laser_)
          UpdateLidar(meas_package);
        break;
      }
      case MeasurementPackage::RADAR: {
        if(use_radar_)
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

  // calculate the square root of lambda + n_aug_
  double coef = sqrt(3.0);  // sqrt(3)
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
    if (fabs(yawd) < 0.0001) {
      Xsig_pred_(0, i) = px + v * cos(yaw) * delta_t
          + 0.5 * dt2 * cos(yaw) * nu_a;
      Xsig_pred_(1, i) = py + v * sin(yaw) * delta_t
          + 0.5 * dt2 * sin(yaw) * nu_a;

    } else {
      Xsig_pred_(0, i) = px + v*(sin(yaw + yawd*delta_t) - sin(yaw))/yawd
          + 0.5 * dt2 * cos(yaw) * nu_a;
      Xsig_pred_(1, i) = py + v*(-cos(yaw + yawd*delta_t) + cos(yaw))/yawd
          + 0.5 * dt2 * sin(yaw) * nu_a;
    }

    Xsig_pred_(2,i) = v + delta_t*nu_a;
    Xsig_pred_(3,i) = yaw + yawd * delta_t + 0.5  * nu_yawdd* dt2;
    Xsig_pred_(4,i) = yawd +  nu_yawdd * delta_t;
  }

  //predict state mean
  x_ = Xsig_pred_ * weights_;

  //predict state covariance matrix
  Xsig_pred_diff_= Xsig_pred_ - x_ * MatrixXd::Ones(1,2*n_aug_ +1);
  //Adjust the angle to between -Pi to Pi
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
    while (Xsig_pred_diff_(3, i) > M_PI)
      Xsig_pred_diff_(3, i) -= 2 * M_PI;
    while (Xsig_pred_diff_(3, i) < -M_PI)
      Xsig_pred_diff_(3, i) += 2 * M_PI;
  }

  //cout <<"Xsig_pred_diff" <<endl<<Xsig_pred_diff_.row(3)<<endl;
  P_ = Xsig_pred_diff_ * weights_.asDiagonal()* Xsig_pred_diff_.transpose();

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

  assert(meas_package.raw_measurements_.size() == 2);
  VectorXd z_pred_mean = x_.head(2);

  MatrixXd S = P_.topLeftCorner(2,2) + R_lidar_;
  MatrixXd Sinv = S.inverse();


  // create Cross-correlation Matrix
  MatrixXd Tc = P_.topLeftCorner(n_x_,2);
  // create Kalman Gain Matrix
  MatrixXd K = Tc * Sinv;
  VectorXd z = meas_package.raw_measurements_;

  x_ += K * (z - z_pred_mean);
  P_ -= K * S * K.transpose();


  //calculate NIS
  double nis_error = (z - z_pred_mean).transpose() * Sinv * (z -z_pred_mean);
  cout << "Lidar NIS error: "<< nis_error<< endl;
  ofs_nis_lidar_<<nis_error<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:Done

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  assert(meas_package.raw_measurements_.size()==3);

  MatrixXd Zsig_pred = MatrixXd(3, 2*n_aug_ +1);
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
    r= max(r,0.001);
    r_dot = (px * v*cos(yaw) + py*v*sin(yaw))/r;
    Zsig_pred.col(i)<< r,phi,r_dot;
  }
  VectorXd z_pred_mean = VectorXd::Zero(3);

  z_pred_mean = Zsig_pred * weights_;

  MatrixXd S = MatrixXd::Zero(3,3);

  //create a Matrix of difference between Zsig to Zsig_mean
  MatrixXd Zdiff = Zsig_pred - z_pred_mean * MatrixXd::Ones(1,2*n_aug_+1);
  //Adjust the angle to between -Pi to Pi
  for(unsigned int i= 0; i< 2*n_aug_ +1;i++){
    while (Zdiff(1, i) > M_PI)
      Zdiff(1, i) -= 2 * M_PI;
    while (Zdiff(1, i) < -M_PI)
      Zdiff(1, i) += 2 * M_PI;
  }

  S = Zdiff * weights_.asDiagonal()*Zdiff.transpose() + R_radar_;

  // create Cross-correlation Matrix
  MatrixXd Tc = Xsig_pred_diff_ * weights_.asDiagonal() * Zdiff.transpose();

  // create Kalman Gain Matrix
  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  VectorXd z = meas_package.raw_measurements_;

  x_ += K*(z - z_pred_mean);
  P_ -= K * S* K.transpose();

  //calculate NIS
  double nis_error = (z - z_pred_mean).transpose() * Sinv * (z -z_pred_mean);
  cout << "Radar NIS error: "<< nis_error<< endl;
  ofs_nis_radar_<<nis_error<<endl;
}
