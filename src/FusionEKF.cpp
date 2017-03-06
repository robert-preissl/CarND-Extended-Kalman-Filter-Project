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

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */

  //set the acceleration noise components
  noise_ax = 5;
  noise_ay = 5;
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

    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    VectorXd x_init = VectorXd(4);
    //set the state with the initial location and zero velocity
		x_init << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    // ekf_.x_ << 1, 1, 1, 1;

    //state covariance matrix P
  	MatrixXd P_init = MatrixXd(4, 4);
  	P_init  << 1, 0, 0, 0,
  		 	       0, 1, 0, 0,
  			       0, 0, 1000, 0,
  			       0, 0, 0, 1000;

    //measurement covariance
    MatrixXd R_init = MatrixXd(2, 2);
    R_init << 0.0225, 0,
              0, 0.0225;

    //measurement matrix
    MatrixXd H_init = MatrixXd(2, 4);
    H_init << 1, 0, 0, 0,
              0, 1, 0, 0;

    //the initial transition matrix F_
    MatrixXd F_init = MatrixXd(4, 4);
    F_init << 1, 0, 1, 0,
              0, 1, 0, 1,
             	0, 0, 1, 0,
             	0, 0, 0, 1;

    //set the process covariance matrix Q to 0 (since dt=0)
    float dt = 0;
  	float dt_2 = dt * dt;
  	float dt_3 = dt_2 * dt;
  	float dt_4 = dt_3 * dt;

    MatrixXd Q_init = MatrixXd(4, 4);
    Q_init <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
            	 0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
            	 dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
            	 0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

    // call the Kalman filter init function to set initial vectors and matrices
    ekf_.Init(x_init, P_init, F_init, H_init, R_init, Q_init);



    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
  } else {
    // Laser updates
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
