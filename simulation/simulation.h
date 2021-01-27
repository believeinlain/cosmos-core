#ifndef SIMQULATION_H
#define SIMQULATION_H

#include <iostream>
#include <vector>
#include "gnuplot_i.hpp"
#include "sph_sim.h"

bool wait_for_key();

class simulation {
private:
	sph_sim SPH;
	param_struct param;
	group_conf_struct group_conf;
	MatrixXd x;
	MatrixXd y;
	MatrixXd u;
	MatrixXd v;
	// Used for plotting the vehicle paths
	vector<double> trackt;
	double plotdt;
	double plott;
	double t0;
	double tf;
	// Gnuplot window
	Gnuplot g1;

	// for ...

	// Loiter circle locations
	// Position [x y]
	MatrixXd lx;
	MatrixXd lR;
	MatrixXd rdx;

	// Make sure to match group_conf.obs_init
	MatrixXd obx;

	void init_simulation();
public:
	simulation();
	void start_simulation();
	void plot_veh(const sph_sim& SPH, const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx);
};



#endif