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
	param_struct param = SPH.get_param();
	group_conf_struct group_conf = SPH.get_group_conf();
	MatrixXd x = SPH.get_x();
	MatrixXd y = SPH.get_y();
	MatrixXd u = SPH.get_u();
	MatrixXd v = SPH.get_v();
	// Used for plotting the vehicle paths
	double trackt = SPH.get_initial_time();
	double plotdt = 0.1;
	double plott = SPH.get_time();
	double t0 = SPH.get_time();
	double tf = 100;
	// for ...

	// Loiter circle locations
	// Position [x y]
	Vector2d lx;
	lx << 28,0;

	// Make sure to match group_conf.obs_init
	MatrixXd obx(5,2);
	obx << 	 7, 0,
			12, 4,
			16, 2,
			 9,-2,
			22,-6;
	
	// Loiter circle radii
	MatrixXd lR;
	if(group_conf.num_loiter > 0) {
		if(SPH.get_time() < 15) {
			// Loiter circle radii
			lR.resize(1,1);
			lR << 5;
		} else {
			// Change the loiter circle radii at each time step - DISABLED
			// lR= ... this is disabled

			// Update the SPH properties
			group_conf.veh_h = 2 * lR.array() * sin(group_conf.num_veh.cast<double>().array() / M_PI);
			SPH.sph_update_properties(param, group_conf);
		}
	} else {
		lR.resize(0,0);
	}

	// Reduced density targets for the vehicles
	// First rd position [x y]
	Matrix2d rdx;
	rdx << 28, 0;

	// Take a SPH timestep
	//SPH.sph_sim_step(rdx,lx,lR);
	
	
	// endfor 
}
