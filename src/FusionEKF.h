#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include <vector>
#include <string>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  int ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Initialize vectors/matrices
  */
  int Initialize(const MeasurementPackage &measurement_pack);

  /**
  * Predict step
  */
  void Predict(const MeasurementPackage &measurement_pack);


  /**
  * Update step
  */
  void Update(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

private:
  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  MatrixXd R_laser_;
  MatrixXd R_radar_;
  MatrixXd H_laser_;
  MatrixXd Hj_;

  //acceleration noise components
  float noise_ax;
  float noise_ay;
};

#endif /* FusionEKF_H_ */
