#include "simulation.h"

simulation::simulation() {
	init_simulation();
}

void simulation::init_simulation() {
	param = SPH.get_param();
	group_conf = SPH.get_group_conf();
	x = SPH.get_x();
	y = SPH.get_y();
	u = SPH.get_u();
	v = SPH.get_v();

	// Used for plotting the vehicle paths
	trackt.push_back(SPH.get_initial_time());
	plotdt = 0.1;	// replot every 100 milliseconds
	plott = SPH.get_time();
	t0 = SPH.get_time();
	tf = 50;	// time final

	lx.resize(1,2);
	rdx.resize(1,2);
	obx.resize(5,2);
}

void simulation::start_simulation() {
	GnuplotPipe gp;
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
				// Change the loiter circle radii at each time step - DISABLED
				// lR= ... this is disabled

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
		x = append_right(x, SPH.get_x());
		y = append_right(y, SPH.get_y());
		u = append_right(u, SPH.get_u());
		v = append_right(v, SPH.get_v());
		trackt.push_back(SPH.get_time());

		if(x.array().isNaN().any()) {
			cout << "Something went wrong, NaN detected in x-positions.";
			throw "Something went wrong, NaN detected in x-positions";
		}

		// Plot
		if(SPH.get_time() >= plott - SPH.get_dt()/10) {
			plot_veh(SPH,x,y,trackt,lx, obx);
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

void simulation::plot_veh(const sph_sim& SPH, const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx) {
	ostringstream o;
	//o << "plot " //<< "(0+word(X,int(t))),(0+word(Y,int(t)))" << ", "			// line ver
	//			 << "(0+word(X,int(t))),(0+word(Y,int(t))) with points";	// point ver
	//gp.sendLine("HI",true);
	
	gp.sendLine(gnuvec(x, "X"), true);
	gp.sendLine(gnuvec(y, "Y"), true);
	gp.sendLine("set parametric", true);
	gp.sendLine("set trange [1:words(X)]; set samples words(X)", true);
	gp.sendLine("unset key", true);			// disable legend
	gp.sendLine("set xrange [-10:40]; set yrange [-10:10]", true);	// set x/y-axis ranges
	//gp.sendLine("plot (0+word(X,int(t))),(0+word(Y,int(t)))", true);
	gp.sendLine("plot (0+word(X,int(t))),(0+word(Y,int(t))) with points", true);
	gp.sendEndOfData();
}



bool wait_for_key () {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return false;
}