#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  cout << "EKF Predict 1 " << endl;
  std::cout << " before predict. x_ = " << x_ << std::endl;
  /**
  TODO:
    * predict the state
  */
printf("MM0 -- x_[0] = %.17g \n", x_(0));
  x_ = F_ * x_;
printf("MM1 -- x_[0] = %.17g / F_[0][0] = %.17g \n", x_(0), F_(0,0));
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  std::cout << " XX3 P_ = " << P_ << " / Q = " << Q_ << std::endl;

  cout << "EKF Predict 2 x_ " << x_ << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  cout << "EKF Update 1 " << endl;
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_; // convert the predicted x state into measurment space
  VectorXd y = z - z_pred;

printf("MM2 -- x_[0] = %.17g / z_pred = %.17g \n", x_(0), z_pred(0));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  std::cout << "AA0 -- z = " << z << " / z_pred = " << z_pred << " / x_ = " << x_ << " / y = " << y << std::endl;
  std::cout << "AA1 -- Ht = " << Ht << " / S = " << S << " / Si = " << Si << " / PHt = " << PHt << " / K = " << K << std::endl;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

printf("MM3 -- x_[0] = %.17g / K[0][0] = %.17g \n", x_(0), K(0,0));
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  cout << "EKF Update 3 " << endl;
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  // Extended differences:
  //  - F_j instead of F
  //  - H_j instead of H

  /*
   You still use Hj_ * x. But Hj_ is only ment to be used for calculating the Kalman Gain.
   To calc the z_predict you need to transform the x_state by your self using h(x).
   */

  Tools tools;

  MatrixXd Hj = tools.CalculateJacobian(x_);
  cout << "Hj = " << Hj << endl;

  cout << "x_ = " << x_ << endl;

  // VectorXd z_pred = H_ * x_; // // convert the predicted x state into measurment space

  double range = sqrt( pow(x_(0),2) + pow(x_(1),2) );
  cout << "range = " << range << endl;
  double bearing = atan(x_(1)/x_(0));
  cout << "bearing = " << range << endl;
  double range_rate =  ((x_(0)*x_(2)+x_(1)*x_(3))/(sqrt( pow(x_(0),2) + pow(x_(1),2) )));
  cout << "range_rate = " << range << endl;
  // MatrixXd zpred(3, 1);
  VectorXd z_pred(3);
  z_pred << range, bearing, range_rate;


  cout << "z_pred = " << z_pred << endl;
  cout << "z = " << z << endl;

  VectorXd y = z - z_pred;
  cout << "y = " << y << endl;

  MatrixXd Hjt = Hj.transpose();
  cout << "Hjt = " << Hjt << endl;

  MatrixXd S = Hj * P_ * Hjt + R_;
  cout << "S = " << S << endl;

  MatrixXd Si = S.inverse();
  cout << "Si = " << Si << endl;

  MatrixXd PHt = P_ * Hjt;
  cout << "PHt = " << PHt << endl;

  MatrixXd K = PHt * Si;
  cout << "K = " << K << endl;

  //new estimate
  x_ = x_ + (K * y);
  cout << "x_ = " << x_ << endl;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
  cout << "P_ = " << P_ << endl;

  cout << "EKF Update 4 " << endl;

}





// done
