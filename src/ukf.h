#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* process noise matrix
  MatrixXd Q_;

  ///* lidar measurement noise matrix
  MatrixXd R_lidar_;
  ///* radar measurement noise matrix
  MatrixXd R_radar_;

  ///*  difference between predicted sigma points to predicted mean, temporal variable
  MatrixXd Xsig_pred_diff_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* file stream to record NIS
  std::ofstream ofs_nis_radar_;
  std::ofstream ofs_nis_lidar_;

public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

private:
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
   * process first measurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessFirstMeasurement(MeasurementPackage meas_package);

  /**
   * unit test to verify the algorithm for augmented sigma points process.
   */
  void test();
  /**
   * generate augmented sigma points
   * @param state current mean of the state vector
   * @param P current Covariant matrix
   * @param Q current process noise matrix
   * @return generated matrix of augmented sigma points
   */
  MatrixXd generateAugmentedSigmaPoints(const VectorXd& state,const MatrixXd& P,const MatrixXd& Q);
  /**
    * predict augmented sigma points
    * @param delta_t elapsed time since last processing
    * @param Xsig_aug Matrix of augmented sigma points before prediction
    * @return Matrix of predicted augmented sigma points
    */
  MatrixXd predictAugmentedSigmaPoints(const double delta_t,const MatrixXd& Xsig_aug);
  /**
    * calculate the mean and Covariant Matrix
    * @param sgima_p Matrix of augmented sigma points
    * @param[out] state mean of augmented sigma points
    * @param[out] Q  process noise matrix
    * @return matrix of deviation of sigma points to mean
    */
  MatrixXd calculatePredictedMeanAndCovariantMatrix(const MatrixXd& sigma_p,VectorXd& state,MatrixXd& P);
};


#endif /* UKF_H */
