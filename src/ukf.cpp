#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <assert.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
//using std::vector;

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
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(VectorXd::Ones(n_x_).asDiagonal());

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;

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

  test();
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
    while(delta_t > 0.1){
      delta_t -= 0.1;
      Prediction(0.1);
    }
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

  // generate augmented sigma points
  MatrixXd Xsig_aug = generateAugmentedSigmaPoints(x_,P_,Q_);
  //predict sigma points

  Xsig_pred_ = predictAugmentedSigmaPoints(delta_t, Xsig_aug);

  //predict state mean
  Xsig_pred_diff_ = calculatePredictedMeanAndCovariantMatrix(Xsig_pred_,x_,P_);
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
    // numerical stability check
    if(fabs(px)<0.005){ // 5mm distance
      px = copysign(0.005,px);
    }
    r = sqrt(px*px + py*py);//r will not be zero
    phi = atan2(py,px); // px,py will not be both zeros
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
  //numerical stability check
  if(fabs(x_[0]) < 0.005)
    x_[0]= copysign(0.005,x_[0]);
}

MatrixXd UKF::generateAugmentedSigmaPoints(const VectorXd& state,const MatrixXd& P,const MatrixXd& Q) {

  // calculate the square root of lambda + n_aug_
  double coef = sqrt(3.0);  // sqrt(3)
  // create augmented mean vector
  VectorXd x_aug_mean = VectorXd::Zero(n_aug_);
  x_aug_mean.head(n_x_) = state;
  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  // create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  P_aug.topLeftCorner(n_x_, n_x_) = P;
  P_aug.bottomRightCorner(2, 2) = Q;
  // create square root matrix of P_aug
  MatrixXd L = P_aug.llt().matrixL();
  // assign each column with x_aug_mean
  Xsig_aug = x_aug_mean * MatrixXd::Ones(1, 2 * n_aug_ + 1);
  Xsig_aug.block(0, 1, n_aug_, n_aug_) += coef * L;
  Xsig_aug.block(0, 1 + n_aug_, n_aug_, n_aug_) -= coef * L;
  return Xsig_aug;
}

MatrixXd UKF::predictAugmentedSigmaPoints(const double delta_t, const MatrixXd& Xsig_aug) {
  //predict sigma points
  double dt2 = delta_t * delta_t;
  // initialize the predict sigma points matrix
  MatrixXd   Xsig_pred = MatrixXd::Zero(n_x_,2 * n_aug_+1);
  for (unsigned int i = 0; i < 2 * n_aug_+1; i++) {
    VectorXd col = Xsig_aug.col(i);
    double px = col[0];
    double py = col[1];
    double v = col[2];
    double yaw = col[3];
    double yawd = col[4];
    double nu_a = col[5];
    double nu_yawdd = col[6];
    if (fabs(yawd) < 0.0001) {
      Xsig_pred(0, i) = px + v * cos(yaw) * delta_t
          + 0.5 * dt2 * cos(yaw) * nu_a;
      Xsig_pred(1, i) = py + v * sin(yaw) * delta_t
          + 0.5 * dt2 * sin(yaw) * nu_a;
    } else {
      Xsig_pred(0, i) = px + v * (sin(yaw + yawd * delta_t) - sin(yaw)) / yawd
          + 0.5 * dt2 * cos(yaw) * nu_a;
      Xsig_pred(1, i) = py + v * (-cos(yaw + yawd * delta_t) + cos(yaw)) / yawd
          + 0.5 * dt2 * sin(yaw) * nu_a;
    }
    Xsig_pred(2, i) = v + delta_t * nu_a;
    Xsig_pred(3, i) = yaw + yawd * delta_t + 0.5 * nu_yawdd * dt2;
    Xsig_pred(4, i) = yawd + nu_yawdd * delta_t;
  }

  return Xsig_pred;
}

MatrixXd UKF::calculatePredictedMeanAndCovariantMatrix(const MatrixXd& XSig,VectorXd& state,MatrixXd& P) {
  //predict state mean
  state = XSig * weights_;
  //create the deviation matrix of predicted sigma points to predicted mean state
  MatrixXd Xsig_diff = XSig - state * MatrixXd::Ones(1, 2 * n_aug_ + 1);
  //Adjust the angle to between -Pi to Pi
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {
    while (Xsig_diff(3, i) > M_PI)
      Xsig_diff(3, i) -= 2 * M_PI;
    while (Xsig_diff(3, i) < -M_PI)
      Xsig_diff(3, i) += 2 * M_PI;
  }
  //get state process Covariance matrix
  P = Xsig_diff * weights_.asDiagonal() * Xsig_diff.transpose();
  return Xsig_diff;
}

void UKF::test()
{
  // test genearte sigma points
  //set example state
   VectorXd x = VectorXd(5);
   x <<   5.7441,
          1.3800,
          2.2049,
          0.5015,
          0.3528;

   //create example covariance matrix
   MatrixXd P = MatrixXd(5, 5);
   P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
           -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
            0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
           -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
           -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

   MatrixXd Q = MatrixXd(2,2);

   Q<<     0.04,       0.0,
           0,          0.04;

   cout << "Test generateAugmentedSigmaPoints...";
   MatrixXd Xsig_aug = generateAugmentedSigmaPoints(x,P,Q);
   MatrixXd Xsig_aug_expected = MatrixXd(7,15);
   Xsig_aug_expected <<
       5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,  5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
         1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,  1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
       2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,  2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
       0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,  0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
       0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528, 0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
            0,        0,        0,        0,        0,        0,  0.34641,        0,        0,        0,        0,        0,        0, -0.34641,        0,
            0,        0,        0,        0,        0,        0,        0,  0.34641,        0,        0,        0,        0,        0,        0, -0.34641;

   assert(Xsig_aug.isApprox(Xsig_aug_expected,0.0001));
   cout << "Passed."<<endl;

   //test predictAugmentedSigmaPoints
   //create example sigma point matrix

   Xsig_aug <<
       5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
         1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
       2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
       0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
       0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
            0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
            0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;


   //create matrix with predicted sigma points as columns

  cout << "Test predictAugmentedSigmaPoints...";
  double delta_t = 0.1; //time diff in sec
  MatrixXd Xsig_pred = predictAugmentedSigmaPoints(delta_t,Xsig_aug);

  MatrixXd Xsig_pred_expected = MatrixXd(5, 15);
  Xsig_pred_expected <<
      5.93553, 6.06251, 5.92217, 5.9415, 5.92361, 5.93516, 5.93705, 5.93553, 5.80832, 5.94481, 5.92935, 5.94553, 5.93589, 5.93401, 5.93553,
      1.48939, 1.44673, 1.66484, 1.49719, 1.508, 1.49001, 1.49022, 1.48939, 1.5308, 1.31287, 1.48182, 1.46967, 1.48876, 1.48855, 1.48939,
      2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.23954, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.17026, 2.2049,
      0.53678, 0.473387, 0.678098, 0.554557, 0.643644, 0.543372, 0.53678, 0.538512, 0.600173, 0.395462, 0.519003, 0.429916, 0.530188, 0.53678, 0.535048,
      0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.387441, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.318159;

  assert(Xsig_pred.isApprox(Xsig_pred_expected,0.0001));
  cout << "Passed."<<endl;


  cout << "Test calculatePredictedMeanAndCovariantMatrix...";

  Xsig_pred <<
           5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
             1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
            2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
           0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
            0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  MatrixXd XSig_diff = calculatePredictedMeanAndCovariantMatrix(Xsig_pred,x,P);
  VectorXd x_expected = VectorXd(5);
  x_expected << 5.93637, 1.49035, 2.20528,  0.536853, 0.353577;
  MatrixXd P_expected = MatrixXd(5,5);
  P_expected <<
      0.00543425, -0.0024053, 0.00341576, -0.00348196, -0.00299378,
      -0.0024053, 0.010845, 0.0014923, 0.00980182, 0.00791091,
      0.00341576, 0.0014923, 0.00580129, 0.000778632, 0.000792973,
      -0.00348196, 0.00980182, 0.000778632, 0.0119238, 0.0112491,
      -0.00299378, 0.00791091, 0.000792973, 0.0112491, 0.0126972;

  assert(x.isApprox(x_expected,0.0001));
  assert(P.isApprox(P_expected,0.0001));
  cout << "Passed."<<endl;
}
