#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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

  // Initial measurement matrix - laser
  H_laser_ << 1,0,0,0,
              0,1,0,0;
  
  // Set the process and measurement noises
  ax_ = 9.0;
  ay_ = 9.0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
      
    // first measurement
    cout << "First measurement" << endl;
    VectorXd x(4);

    previous_timestamp_ = measurement_pack.timestamp_;
    
    // Covariance matrix
    MatrixXd P(4,4);
    P << 1,0,0,0,
         0,1,0,0,
         0,0,1000,0,
         0,0,0,1000;
    
    // Transition matrix
    MatrixXd F(4,4);
    F << 1,0,0,0,
         0,1,0,0,
         0,0,1,0,
         0,0,0,1;
    
    MatrixXd Q(4,4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      x << rho*cos(phi), rho*sin(phi), 0.0f, 0.0f;
      Hj_ = tools.CalculateJacobian(x);
      ekf_.Init(x, P, F, Hj_, R_radar_, Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0f, 0.0f;
      ekf_.Init(x, P, F, H_laser_, R_laser_,Q);
    }
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1,0,dt,0,
              0,1,0,dt,
              0,0,1,0,
              0,0,0,1;
  
  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double dt4 = dt3*dt;
  double dt4_4 = dt4/4.0;
  double dt3_2 = dt3/2.0;
  ekf_.Q_ << dt4_4*ax_,0,dt3_2*ax_,0,
               0,dt4_4*ay_,0,dt3_2*ay_,
               dt3_2*ax_,0,dt2*ax_,0,
               0,dt3_2*ay_,0,dt2*ay_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
