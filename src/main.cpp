
// http://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/
// http://math.stackexchange.com/questions/840662/an-explanation-of-the-kalman-filter
//
//  PREDICT STEP:
//
//   1. compute next state:
//    use F .. prediction matrix, which gives us our next state:  x_k = F_k * x_k-1  (x = (p,v))
//
//   2. update covariance matrix
//    P_k = F_k * P_k−1 * FT_k
//
//   3. model the uncertainty associated with the “world” (i.e. things we aren’t keeping track of)
//    by adding some new uncertainty after every prediction step:
//   then, we produce a new Gaussian blob, with a different covariance (but the same mean):
//    by simply adding Qk, giving our complete expression for the prediction step:
//
//   -> we add Q_k to P_k (and add B_k*u_k to x_k)
//
//    -> P_k = F_k * P_k−1 * FT_k + Q_k
//
//
//   We have a fuzzy estimate of where our system might be, given by x̂ k and Pk.
//
//   !! What happens when we get some data from our sensors?
//
//
//  Each sensor tells us something indirect about the state— in other words,
//   the sensors operate on a state and produce a set of readings.
//
//  Notice that the units and scale of the reading might not be the same as the units
//   and scale of the state we’re keeping track of.
//  ==> We’ll model the sensors with a matrix, H_k.
//
//
//  From each reading we observe, we might guess that our system was in a particular state.
//   But because there is uncertainty, some states are more likely than others to have have
//   produced the reading we saw.
//
//  terms:
//     R_k .. covariance of this uncertainty (i.e. of the sensor noise)
//     z_k .. the reading we observed, or the mean of the distribution, the k'th measurement
//
//  A) -> So now we have two Gaussian blobs (in measurement space):
//       - One surrounding the mean of our transformed prediction (E = H_k * P_k * HT_k);  H maps from state space to measurement space
//       - one surrounding the actual sensor reading we got.
//
//  b) -> And, we have two estimates (in measurement space): 1) H_k * x_k  + 2) z_k (the sensor reading)
//
//
// So what’s our new most likely state?
//
//  !! If we have two probabilities and we want to know the chance that both are true, we just multiply them together.
//
//  What we’re left with is the overlap, the region where both blobs are bright/likely. And it’s a lot more precise than either of our previous estimates. The mean of this distribution is the configuration for which both estimates are most likely, and is therefore the best guess of the true configuration given all the information we have.
//  As it turns out, when you multiply two Gaussian blobs with separate means and covariance matrices, you get a new Gaussian blob with its own mean and covariance matrix!
//
//
//   Putting it all together for the UPDATE STEP:
//
//    We have two distributions: (in measurement space)
//     1.) The predicted measurement : (μ0,Σ0) = (H_k * x̂_k, H_k * P_k * HT_k)
//     2.) the observed measurement  : (μ1,Σ1) = (z_k, R_k)
//
//
//   K′ = P_k * HT_k * (H_k * P_k * HT_k + R_k)^-1   .. Kalman gain
//   P_k' = P_k - K′ * H_k * P_k
//   x_k' = x_k + K′(z_k - H_k * x̂_k)
//
//
//   And that’s it! x̂ ′k is our new best estimate, and we can go on and feed it (along with P′k ) back into another round of predict or update as many times as we like.


#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "FusionEKF.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {

    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long timestamp;

    // reads first element from the current line
    iss >> sensor_type;
    if (sensor_type.compare("L") == 0) {
      // LASER MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float x;
      float y;
      iss >> x;
      iss >> y;
      meas_package.raw_measurements_ << x, y;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // RADAR MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float theta;
      float ro_dot;
      iss >> ro;
      iss >> theta;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, theta, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  // Create a Fusion EKF instance
  FusionEKF fusionEKF;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  //Call the EKF-based fusion
  size_t N = measurement_pack_list.size();
  for (size_t k = 0; k < N; ++k) {
    // start filtering from the second frame (the speed is unknown in the first frame)
    if(fusionEKF.ProcessMeasurement(measurement_pack_list[k]) < 0) {
      continue; // if we got a negative responese code (e.g., 0 values), we try another measurment to initialize (only return <0 in init phase)
    }
    // output the estimation
    out_file_ << fusionEKF.ekf_.x_(0) << "\t";
    out_file_ << fusionEKF.ekf_.x_(1) << "\t";
    out_file_ << fusionEKF.ekf_.x_(2) << "\t";
    out_file_ << fusionEKF.ekf_.x_(3) << "\t";

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // output the estimation
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // p1_meas
      out_file_ << ro * sin(phi) << "\t"; // ps_meas
    }

    // output the ground truth packages
    out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(3) << "\n";

    estimations.push_back(fusionEKF.ekf_.x_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }

  // compute the accuracy (RMSE)
  Tools tools;
  cout << "Accuracy - RMSE:" << endl << tools.CalculateRMSE(estimations, ground_truth) << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  return 0;
}
