#include <iostream>
#include <numeric>
#include <vector>
#include "Eigen/Dense"
#include "sph_sim.h"
 
using namespace Eigen;
using namespace std;

int main()
{
	sph_sim SPH;
	MatrixXd x = SPH.get_x();
	MatrixXd y = SPH.get_y();
	MatrixXd u = SPH.get_u();
	MatrixXd v = SPH.get_v();
	// Used for plotting the vehicle paths
	double trackt = SPH.get_initial_time();
}
