#ifndef SIMQULATION_H
#define SIMQULATION_H

#include <iostream>
#include <chrono>         // std::chrono::seconds
#include <thread>         // std::this_thread::sleep_for
#include <vector>
#include "gnuplot.h"
#include "sph.h"

#include "agent/agentclass.h"


/// Agents will send world controller their sph-updated states
int32_t send_world_new_state(string &request, string &response, Agent *agent);

class simulation {
private:
	/// Pointer to simulation control agent
	Agent *agent;
	/// SPH object
	sph_sim SPH;
	/// Simulation parameters
	param_struct param;
	/// Simulation group configurations
	group_conf_struct group_conf;
	/// History of x positions of particles
	Eigen::MatrixXd x;
	/// History of y positions of particles
	Eigen::MatrixXd y;
	/// History of u velocities of particles
	Eigen::MatrixXd u;
	/// History of v velocities of particles
	Eigen::MatrixXd v;
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
	Eigen::MatrixXd lx;
	/// Loiter circle radii
	Eigen::MatrixXd lR;
	/// Reduced density target locations [x y]
	Eigen::MatrixXd rdx;

	/// Make sure to match group_conf.obs_init
	Eigen::MatrixXd obx;

	/// Past number of time steps to track positions of
	int trackMax;

	/// Track the index of the head of the matrix
	int thead;

	/// Number of agents in the simulation
	int num_agents;

	/// Initialize the simulation parameters
	void init_simulation(int tf=100);



	/// Create a visualization for the simulation using gnuplot
	/**
	@param	x			Matrix of history of x positions of all SPH particles
	@param	y			Matrix of history of y positions of all SPH particles
	@param	trackt		Vector of time steps
	@param	lx			Matrix containing [x y] positions of the loiter circles
	@param	obx			Matrix containing [x y] positions of the obstacles
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	void plot_veh(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y, const vector<double>& trackt, const Eigen::MatrixXd& lx, const Eigen::MatrixXd& obx);

	/// Plot the SPH particles
	/**
	@param	x			Matrix of history of x positions of all SPH particles
	@param	y			Matrix of history of y positions of all SPH particles
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	void plot_points(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y);

	/// Plot the loiter circles
	/**
	@param	lx			Matrix containing [x y] positions of the loiter circles
	@return n/a
	*/
	void plot_lx(const Eigen::MatrixXd& lx);

	/// Display trail for each particle
	/**
	@param	x			Matrix of history of x positions of all SPH particles
	@param	y			Matrix of history of y positions of all SPH particles
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	void plot_trails(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y);

	/// Convert history of particle positions into a gnuplot-parsable line
	/**
	@param	pos			Matrix of history of positions of SPH particle
	@param	varname		Variable name to assign the line. Used by gnuplot.
	@param	thead		Index position of the head of the vectors
	@return n/a
	*/
	string gnutrail(const Eigen::RowVectorXd& pos, const string& varname);
	


public:
	/// Constructor
	simulation(Agent *agent, int tf=100);
	/// Start the simulation loop
	void start_simulation();
	/// Check that all agents in the simulation are running
	bool all_sim_agents_running();
	/// Initialize sim agents
	void init_sim_agents(bool test=false);
	/// Send a request to all agents in the simulation
	vector<string> send_req_to_all_agents(string request);
	



	/// Get initial time
	double get_initial_time();
};



#endif