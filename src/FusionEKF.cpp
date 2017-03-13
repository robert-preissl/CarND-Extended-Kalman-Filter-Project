#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXf;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0;

//ekf_.x_ = VectorXd(4); // TODO: otherwise NPE
//  ekf_.x_ << 0, 0, 0, 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4); // for radar?

  // covariance of the measurement uncertainty (i.e. of the sensor noise)
  R_radar_ << 0.0225, 0, 0,
              0, 0.0225, 0,
              0, 0, 0.0225;

  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //set the acceleration noise components
  noise_ax = 15;
  noise_ay = 15;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

int FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
   cout << " ZZ0 -- measurement_pack.raw_measurements_ = " << measurement_pack.raw_measurements_ << endl;

// TODO: move into sep function
  if (!is_initialized_) {

    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF initialization " << endl;
    VectorXd x_init = VectorXd(4);

    // state covariance matrix P
  	MatrixXd P_init = MatrixXd(4, 4);
  	P_init  << 1, 0, 0, 0,
  		 	       0, 1, 0, 0,
  			       0, 0, 1000, 0, // does 1000 mean that we are super certain of the initial velocity??
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
    cout << " ZZ1 -- measurement_pack.raw_measurements_ = " << measurement_pack.raw_measurements_ << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates to set init state.
      */
      VectorXd position_polar_coords(2);
      VectorXd velocity_polar_coords(2);
      VectorXd position_cartesian_coords(2);
      VectorXd velocity_cartesian_coords(2);

      position_polar_coords << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
      velocity_polar_coords << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[2];

      Tools tools;
      position_cartesian_coords = tools.ConvertPolarToCartesian(position_polar_coords);
      velocity_cartesian_coords = tools.ConvertPolarToCartesian(velocity_polar_coords);

      // TODO: take out
    /*  float x_cart = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
      float y_cart = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);

      if (x_cart == 0 || y_cart == 0){
        return;
      }
      x_init << x_cart, y_cart, 0, 0;*/


      // TODO: init vel to 0
      x_init << position_cartesian_coords(0), position_cartesian_coords(1), 0, 0;//  velocity_cartesian_coords(1), velocity_cartesian_coords(1);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "A. measurement_pack.raw_measurements_ = " << measurement_pack.raw_measurements_ << endl;
      if (measurement_pack.raw_measurements_(0) == 0 || measurement_pack.raw_measurements_(1) == 0){
        cout << "A2" << endl;
        return -1;
      }
      cout << "A1" << endl;
      //set the state with the initial location and zero velocity
      x_init << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    printf("SS0 -- x_init[0] = %.17g \n", x_init(0));
    cout << "x_init = " << x_init << endl;

    if (x_init(0) == 0 || x_init(1) == 0) {
      return -1;
    }

    // call the Kalman filter init function to set initial vectors and matrices
    ekf_.Init(x_init, P_init, F_init, H_init, R_init, Q_init);

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    cout << "EKF initialization done " << endl;

    return 0;
  } // end Initialization

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

// TODO: move in function. same for update

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

   //compute the time elapsed between the current and previous measurements
   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
   cout << " dt = " << dt << endl;
   printf("DT0 -- dt = %.17g \n", dt);
   previous_timestamp_ = measurement_pack.timestamp_;

   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;

   // update the state transition F matrix so that the time is integrated
   ekf_.F_(0, 2) = dt; // to get p_x = p_x(k-1) + dt * v_x
   ekf_.F_(1, 3) = dt;

   // update the process covariance matrix Q
   ekf_.Q_ = MatrixXd(4, 4);
   ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
printf("DT1 -- dt = %.17g / noise_ax =  %.17g \n", dt, noise_ax);


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
    // !! convert from polar to cartesian
    cout << " radar measurement: " << measurement_pack.raw_measurements_ << endl;
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates:
    cout << " laser measurement: " << measurement_pack.raw_measurements_ << endl;
    //ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  printf("RR0 -- x_[0] = %.17g \n", ekf_.x_(0));
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

  return 0;
}






// done
