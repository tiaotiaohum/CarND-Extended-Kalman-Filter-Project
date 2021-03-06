#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"


// no actual Kalman filter equations here
// initializing variables, initializing the Kalman filters
// prepare the Q and F matrices for the prediction step
// then calling functions that implement the prediction step or update step

// initialize variables and matrices (x, F, H_laser, H_jacobian, P, etc.)
// initialize the Kalman filter position vector with the first sensor measurements
// modify the F and Q matrices prior to the prediction step based on the elapsed time between measurements
// call the update step for either the lidar or radar sensor measurement. Because the update step for lidar and radar are slightly different, there are different functions for updating lidar and radar.


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  // measurement matrix for laser (to convert state vector to measurement dimension)
  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;
  
  // state/object covarience matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
  
  // state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;
  
  // create empty process covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);
  

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}
// 
// Initialization of the Kalman filter, doing the prediction and update steps of the Kalman filter.
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float rho = measurement_pack.raw_measurements_(0);
      float theta = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);
      
      ekf_.x_(0) = rho * cos(theta);
      ekf_.x_(1) = rho * sin(theta);
      ekf_.x_(2) = rho_dot * cos(theta);
      ekf_.x_(3) = rho_dot * sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }

    previous_timestamp_ = measurement_pack.timestamp_;
//     previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
// The ekf_ variable is an instance of the KalmanFilter class. We will use ekf_ to store  Kalman filter variables (x, P, F, H, R, Q) and call the predict and update functions.
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - is in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  const float noise_ax = 9;
  const float noise_ay = 9;
  
  // Update process noise covariance matrix
  ekf_.Q_ <<  dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
              0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
              dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
              0, dt_3/2 * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  //convert the time elapsed between the current and previous measurements
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    MatrixXd Hj_temp(3,4);
    Hj_temp= Tools().CalculateJacobian(ekf_.x_);
    if(!Hj_temp.isZero()){
      ekf_.H_ = Hj_temp;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    
  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}