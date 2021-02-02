#ifndef SIMQULATION_H
#define SIMQULATION_H

#include <iostream>
#include <chrono>         // std::chrono::seconds
#include <thread>         // std::this_thread::sleep_for
#include <vector>
#include "gnuplot.h"
#include "sph.h"


class simulation {
private:
	/// SPH object
	sph_sim SPH;
	/// Simulation parameters
	param_struct param;
	/// Simulation group configurations
	group_conf_struct group_conf;
	/// History of x positions of particles
	MatrixXd x;
	/// History of y positions of particles
	MatrixXd y;
	/// History of u velocities of particles
	MatrixXd u;
	/// History of v velocities of particles
	MatrixXd v;
	/// Tracked time steps
	vector<double> trackt;
	/// Timestep for plot
	double plotdt;
	/// Plot time
	double plott;
	/// Initial time
	double t0;
	/// Final time, when the simulation ends
	double tf;

	/// Gnuplot window
	GnuplotPipe gp;

	/// Loiter circle locations [x y]
	MatrixXd lx;
	/// Loiter circle radii
	MatrixXd lR;
	/// Reduced density target locations [x y]
	MatrixXd rdx;

	/// Make sure to match group_conf.obs_init
	MatrixXd obx;

	/// Past number of time steps to track positions of
	int trackMax;

	/// Initialize the simulation parameters
	void init_simulation();

	/// Create a visualization for the simulation using gnuplot
	/**
	@param	x			Matrix of history of x positions of all SPH particles
	@param	y			Matrix of history of x positions of all SPH particles
	@param	trackt		Vector of time steps
	@param	lx			Matrix containing [x y] positions of the loiter circles
	@param	obx			Matrix containing [x y] positions of the obstacles
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	void plot_veh(const MatrixXd& x, const MatrixXd& y, const vector<double>& trackt, const MatrixXd& lx, const MatrixXd& obx, const int& thead);

	/// Plot the SPH particles
	/**
	@param	x			Matrix of history of x positions of all SPH particles
	@param	y			Matrix of history of x positions of all SPH particles
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	void plot_points(const MatrixXd& x, const MatrixXd& y, const int& thead);

	/// Plot the loiter circles
	/**
	@param	lx			Matrix containing [x y] positions of the loiter circles
	@return n/a
	*/
	void plot_lx(const MatrixXd& lx);

	/// Display trail for each particle
	/**
	@param	x			Matrix of history of x positions of all SPH particles
	@param	y			Matrix of history of x positions of all SPH particles
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	void plot_trails(const MatrixXd& x, const MatrixXd& y, const int& thead);

	/// Convert history of particle positions into a gnuplot-parsable line
	/**
	@param	pos			Matrix of history of positions of SPH particle
	@param	varname		Variable name to assign the line. Used by gnuplot.
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	string gnutrail(const RowVectorXd& pos, const string& varname, const int& thead);
	


public:
	/// Constructor
	simulation();
	/// Start the simulation loop
	void start_simulation();
};



#endif