#include "simulation.h"

simulation::simulation() {
	init_simulation();
}

void simulation::init_simulation() {
	param = SPH.get_param();
	group_conf = SPH.get_group_conf();
	// n x 100 matrix
	trackMax = 100;
	x = SPH.get_x() * MatrixXd::Ones(1,trackMax);
	y = SPH.get_y() * MatrixXd::Ones(1,trackMax);
	u = SPH.get_u() * MatrixXd::Ones(1,trackMax);
	v = SPH.get_v() * MatrixXd::Ones(1,trackMax);

	// Used for plotting the vehicle paths
	trackt = vector<double>(100);
	trackt[0] = SPH.get_initial_time();
	plotdt = 0.1;	// replot every 100 milliseconds
	plott = SPH.get_time();
	t0 = SPH.get_time();
	tf = 100;	// time final, when to stop simulation

	lx.resize(1,2);
	lx << 28,0;
	rdx.resize(1,2);
	// Make sure to match group_conf.obs_init
	obx.resize(5,2);
	obx << 	 7, 0,
			12, 4,
			16, 2,
			 9,-2,
			22,-6;
}

void simulation::start_simulation() {
	GnuplotPipe gp;
	int thead = 0; // where in the tracking matrix the head is

	for(double t = t0; t < tf; t += SPH.get_dt()) {
		// Loiter circle locations
		// Position [x y]
		// Slow down the simulation a bit
		this_thread::sleep_for (chrono::milliseconds(4));

		
		
		
		// Loiter circle radii
		if(group_conf.num_loiter > 0) {
			// arbitrary value 15 chosen for time value, when to update SPH properties
			if(SPH.get_time() < 15) {
				// Loiter circle radii
				lR.resize(1,1);
				lR << 5;
			} else if(SPH.get_time() < 25) {
				// Change the loiter circle radii after 15 seconds
				lR.resize(1,1);
				lR << 2;

				// Update the SPH properties
				group_conf.veh_h = 2 * lR.array() * sin(group_conf.num_veh.cast<double>().array() / M_PI);
				SPH.sph_update_properties(param, group_conf);
			} else if(SPH.get_time() < 35) {
				// Change the loiter circle radii after 25 seconds
				lR.resize(1,1);
				lR << 9.5;

				// Update the SPH properties
				group_conf.veh_h = 2 * lR.array() * sin(group_conf.num_veh.cast<double>().array() / M_PI);
				SPH.sph_update_properties(param, group_conf);
			} else if(SPH.get_time() < 90) {
				// Change the loiter circle location after 35 seconds
				lx(0) = lx(0) - 0.005;

				// Update the SPH properties
				group_conf.veh_h = 2 * lR.array() * sin(group_conf.num_veh.cast<double>().array() / M_PI);
				SPH.sph_update_properties(param, group_conf);
			} else {

			}
		} else {
			lR.resize(0,0);
		}

		// Reduced density targets for the vehicles
		// First rd position [x y]
		rdx << 28, 0; // NOTE: this doesn't need to change?

		// Take a SPH timestep
		SPH.sph_sim_step(rdx,lx,lR);

		// Keep track of vehicle paths (or at least the last 100 points)
		x.col(thead) = SPH.get_x();
		y.col(thead) = SPH.get_y();
		u.col(thead) = SPH.get_u();
		v.col(thead) = SPH.get_v();
		trackt[thead] = SPH.get_time();

		if(x.array().isNaN().any()) {
			cout << "Something went wrong, NaN detected in x-positions.";
			throw "Something went wrong, NaN detected in x-positions";
		}

		// Plot
		if(SPH.get_time() >= plott - SPH.get_dt()/10) {
			plot_veh(x,y,trackt,lx, obx, thead);
			plott = plott + plotdt;
		}

		thead = (thead + 1) % trackMax;
	}
}

string gnuvec(const MatrixXd& mat, const string& varname) {
	ostringstream o;
	o << varname << "=\"";
	for(auto it : mat(all,last))
		o << it << " ";
	o << "\"";
	return o.str();
}

void simulation::plot_veh(const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx, const int& thead) {
	gp.sendLine("reset", true);
	gp.sendLine("set title \"Smoothed Particle Hydrodynamics for Agent Control\\n{/*0.85Time = " +to_string(trackt[thead]) + "}\" font \"Arial,16\"", true);
	gp.sendLine("set parametric", true);
	plot_points(x,y,thead);
	plot_lx(lx);
	plot_trails(x,y,thead);
}

void simulation::plot_points(const MatrixXd& x, const MatrixXd& y, const int& thead) {
	for(int i = 0; i < x.rows(); ++i) {
		// Vehicle
		if(i < SPH.get_nveh()) {
			// Vehicles are empty circles
			gp.sendLine("set label at " + to_string(x(i,thead)) + "," + to_string(y(i,thead)) + " point pointtype 6 pointsize 1 lt rgb \"royalblue\"", true);
		}
		// Obstacle
		else if(i < SPH.get_nveh() + SPH.get_nobs()) {
			// Obstacles are empty squares
			gp.sendLine("set label at " + to_string(x(i,thead)) + "," + to_string(y(i,thead)) + " point pointtype 4 pointsize 2 lt rgb \"salmon\"", true);
		}
		// Reduced density particle
		else {
			// Reduced density particles are stars
			gp.sendLine("set label at " + to_string(x(i,thead)) + "," + to_string(y(i,thead)) + " point pointtype 3 pointsize 2 lt rgb \"goldenrod\"", true);
		}
	}
}

// Plot the loiter circle
void simulation::plot_lx(const MatrixXd& lx) {
	if(lx.size() != 0) {
		for(auto row : lx.rowwise()) {
			gp.sendLine("set label at " + to_string(row(0)) + "," + to_string(row(1)) + " point pointtype 3 pointsize 2 lt rgb \"goldenrod\"", true);
		}
	}
}

// Display trail for each particle
void simulation::plot_trails(const MatrixXd& x, const MatrixXd& y, const int& thead) {
	for(int i = 0; i < x.rows(); ++i) {
		gp.sendLine(gnutrail(x.row(i),"X"+to_string(i),thead), true);
		gp.sendLine(gnutrail(y.row(i),"Y"+to_string(i),thead), true);
	}
	gp.sendLine("set trange [1:words(X0)]; set samples words(X0)", true);
	gp.sendLine("unset key", true);
	gp.sendLine("set xrange [-10:40]; set yrange [-10:10]; set size ratio -1", true);
	ostringstream o;
	o << "plot ";
	for(int i = 0; i < x.rows(); ++i) {
		o << "(0+word(X" << to_string(i) << ",int(t))),(0+word(Y" << to_string(i) << ",int(t))) lt rgb \"royalblue\"";
		if(i+1 < x.rows()) {
			o << ", ";
		}
	}
	gp.sendLine(o.str(), true);
	gp.sendEndOfData();
}

// Display trail for each particle
string simulation::gnutrail(const RowVectorXd& pos, const string& varname, const int& thead) {
	ostringstream o;
	o << varname << "=\"";
	for(int i = 0; i < pos.size(); ++i) {
		int tidx = (thead+i+1) % trackMax;
		o << pos(tidx) << " ";
	}
	o << "\"";
	return o.str();
}
