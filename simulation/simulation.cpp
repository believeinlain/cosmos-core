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
	trackt.push_back(SPH.get_initial_time());
	plotdt = 0.1;	// replot every 100 milliseconds
	plott = SPH.get_time();
	t0 = SPH.get_time();
	tf = 30;	// time final, when to stop simulation

	lx.resize(1,2);
	rdx.resize(1,2);
	obx.resize(5,2);
}

void simulation::start_simulation() {
	GnuplotPipe gp;
	int thead = 0; // where in the tracking matrix the head is

	for(double t = t0; t < tf; t += SPH.get_dt()) {
		// Loiter circle locations
		// Position [x y]
		lx << 28,0;
		// Slow down the simulation a bit
		this_thread::sleep_for (chrono::milliseconds(5));

		// Make sure to match group_conf.obs_init
		obx << 	 7, 0,
				12, 4,
				16, 2,
				 9,-2,
				22,-6;
		
		// Loiter circle radii
		if(group_conf.num_loiter > 0) {
			// arbitrary value 15 chosen for time value, when to update SPH properties
			if(SPH.get_time() < 15) {
				// Loiter circle radii
				lR.resize(1,1);
				lR << 5;
			} else {
				// Change the loiter circle radii after 15 seconds
				lR.resize(1,1);
				lR << 2;

				// Update the SPH properties
				group_conf.veh_h = 2 * lR.array() * sin(group_conf.num_veh.cast<double>().array() / M_PI);
				SPH.sph_update_properties(param, group_conf); // This gets called on the crash
			}
		} else {
			lR.resize(0,0);
		}

		// Reduced density targets for the vehicles
		// First rd position [x y]
		rdx << 28, 0;

		// Take a SPH timestep
		SPH.sph_sim_step(rdx,lx,lR);

		// Keep track of vehicle paths (or at least the last 100 points)
		x.col(thead) = SPH.get_x();
		y.col(thead) = SPH.get_y();
		u.col(thead) = SPH.get_u();
		v.col(thead) = SPH.get_v();
		trackt.push_back(SPH.get_time());
		thead = (thead + 1) % trackMax;

		if(x.array().isNaN().any()) {
			cout << "Something went wrong, NaN detected in x-positions.";
			throw "Something went wrong, NaN detected in x-positions";
		}

		// Plot
		if(SPH.get_time() >= plott - SPH.get_dt()/10) {
			plot_veh(SPH,x,y,trackt,lx, obx, thead);
			plott = plott + plotdt;
		}

	}
	//wait_for_key();

}

string gnuvec(const MatrixXd& mat, const string& varname) {
	ostringstream o;
	o << varname << "=\"";
	for(auto it : mat(all,last))
		o << it << " ";
	o << "\"";
	return o.str();
}

void simulation::plot_veh(const sph_sim& SPH, const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx, const int& thead) {
	ostringstream o;
	//o << "plot " //<< "(0+word(X,int(t))),(0+word(Y,int(t)))" << ", "			// line ver
	//			 << "(0+word(X,int(t))),(0+word(Y,int(t))) with points";	// point ver
	//gp.sendLine("HI",true);
	
	// Set variables "X" and "Y" to be the xy positions of the vehicles
	/*gp.sendLine(gnuvec(x, "X"), true);
	gp.sendLine(gnuvec(y, "Y"), true);
	// Plot mode set to 2D
	gp.sendLine("set parametric", true);
	// Configure variable t to iterate through the xy positions
	gp.sendLine("set trange [1:words(X)]; set samples words(X)", true);
	// Disable legend in the corner
	gp.sendLine("unset key", true);
	// Set x/y-axis ranges
	gp.sendLine("set xrange [-10:40]; set yrange [-10:10]", true);


	//gp.sendLine("plot (0+word(X,int(t))),(0+word(Y,int(t)))", true);
	gp.sendLine("plot (0+word(X,int(t))),(0+word(Y,int(t))) with points linetype 6", true);
	gp.sendEndOfData();*/
	plot_trails(x,y,thead);
}

string simulation::gnutrail(const RowVectorXd& pos, const string& varname, const int& thead) {
	ostringstream o;
	o << varname << "=\"";
	for(int i = 0; i < pos.size(); ++i) {
		int tidx = (thead+i) % trackMax;
		o << pos(tidx) << " ";
	}
	o << "\"";
	return o.str();
}

void simulation::plot_trails(const MatrixXd& x, const MatrixXd& y, const int& thead) {
	for(int i = 0; i < x.rows(); ++i) {
		gp.sendLine(gnutrail(x.row(i),"X"+to_string(i),thead), true);
		gp.sendLine(gnutrail(y.row(i),"Y"+to_string(i),thead), true);
	}

	gp.sendLine("set parametric", true);
	gp.sendLine("set trange [1:words(X0)]; set samples words(X0)", true);
	gp.sendLine("unset key", true);
	gp.sendLine("set xrange [-10:40]; set yrange [-10:10]", true);
	ostringstream o;
	o << "plot ";
	for(int i = 0; i < x.rows(); ++i) {
		o << "(0+word(X" << to_string(i) << ",int(t))),(0+word(Y" << to_string(i) << ",int(t)))";
		if(i+1 < x.rows()) {
			o << ", ";
		}
	}
	gp.sendLine(o.str(), true);
	gp.sendEndOfData();
}