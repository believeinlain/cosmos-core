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
	vector<MatrixXd> x { SPH.get_x() };
	vector<MatrixXd> y { SPH.get_y() };
	vector<MatrixXd> u { SPH.get_u() };
	vector<MatrixXd> v { SPH.get_v() };
	// Used for plotting the vehicle paths
	vector<double> trackt { SPH.get_initial_time() };
	double plotdt = 0.1;
	double plott = SPH.get_time();
	double t0 = SPH.get_time();
	double tf = 100;
	// for ...

	// Loiter circle locations
	// Position [x y]
	MatrixXd lx(1,2);
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
	MatrixXd rdx(1,2);
	rdx << 28, 0;

	// Take a SPH timestep
	SPH.sph_sim_step(rdx,lx,lR);

	// Keep track of vehicle paths (or at least the last 100 points)
	x.push_back(SPH.get_x());
	y.push_back(SPH.get_y());
	u.push_back(SPH.get_u());
	v.push_back(SPH.get_v());
	trackt.push_back(SPH.get_time());
	
	for(auto it : x) {
		cout << it << endl << endl;
	}
	for(auto it : y) {
		cout << it << endl << endl;
	}
	for(auto it : u) {
		cout << it << endl << endl;
	}
	for(auto it : v) {
		cout << it << endl << endl;
	}
	for(auto it : trackt) {
		cout << it << endl << endl;
	}
	
	
	// endfor 
}
