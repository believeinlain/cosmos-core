#include <iostream>
#include <numeric>
#include <vector>
#include "Eigen/Dense"
#include "sph_sim.h"
#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communikation with Gnuplot
 
using namespace Eigen;
using namespace std;

void plot_veh(const sph_sim& SPH, const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx);
void wait_for_key(); 	// Programm halts until keypress
Gnuplot g1;

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
	vector<double> trackt { SPH.get_initial_time() };
	double plotdt = 0.1;
	double plott = SPH.get_time();
	double t0 = SPH.get_time();
	double tf = 100;
	g1.showonscreen(); // window output

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
	x = append_right(x, SPH.get_x());
	y = append_right(y, SPH.get_y());
	u = append_right(u, SPH.get_u());
	v = append_right(v, SPH.get_v());
	trackt.push_back(SPH.get_time());
	
	/*cout << "x:\n" << x << endl << endl;
	cout << "y:\n" << y << endl << endl;
	cout << "u:\n" << u << endl << endl;
	cout << "v:\n" << v << endl << endl;
	for(auto it : trackt) {
		cout << it << " ";
	} cout << endl << endl;*/

	if(x.array().isNaN().any()) {
		cout << "Something went wrong, NaN detected in x-positions.";
		throw "Something went wrong, NaN detected in x-positions";
	}

	// Plot
	if(SPH.get_time() >= plott - SPH.get_dt()/10) {
		plot_veh(SPH,x,y,trackt,lx, obx);
		plott = plott + plotdt;
	}
	
	// endfor 
}

void plot_veh(const sph_sim& SPH, const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx) {
	VectorXd v1 = x(all,last);
	vector<double> vx;
	vx.resize(v1.size());
	VectorXd::Map(&vx[0], v1.size()) = v1;

	v1 = y(all,last);
	vector<double> vy;
	vy.resize(v1.size());
	VectorXd::Map(&vy[0], v1.size()) = v1;

	g1.reset_plot();
	cout << endl << endl << "*** user-defined lists of points (x,y)" << endl;
	g1.set_grid();
	g1.set_style("points").plot_xy(vx,vy,"user-defined points 2d");
}



void wait_for_key ()
{
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
    return;
}