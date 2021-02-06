#include "simulation.h"


/// Constructor
simulation::simulation(Agent *agent, int tf/*=100*/) : agent(agent), tf(tf) {
	init_simulation();
}

/// Initialize the simulation parameters
void simulation::init_simulation(int tf/*=100*/) {
	param = SPH.get_param();
	group_conf = SPH.get_group_conf();
	// n x 100 matrix
	trackMax = 100;
	thead = 0;
	x = SPH.get_x() * Eigen::MatrixXd::Ones(1,trackMax);
	y = SPH.get_y() * Eigen::MatrixXd::Ones(1,trackMax);
	u = SPH.get_u() * Eigen::MatrixXd::Ones(1,trackMax);
	v = SPH.get_v() * Eigen::MatrixXd::Ones(1,trackMax);

	// Used for plotting the vehicle paths
	trackt = vector<double>(100);
	trackt[0] = SPH.get_initial_time();
	plotdt = 0.1;	// replot every 100 milliseconds
	plott = SPH.get_time();
	t0 = SPH.get_time();

	// Loiter circle position
	lx.resize(1,2);
	lx << 28,0;
	// Attractors
	rdx.resize(1,2);
	// Obstacle positions
	// Make sure to match group_conf.obs_init
	obx.resize(5,2);
	obx << 	 7, 0,
			12, 4,
			16, 2,
			 9,-2,
			22,-6;
	
	// Initialize agents
	num_agents = 9;
}
string request = "are_you_out_there";
string response = "";
int32_t are_you_out_there(string & request, string &response, Agent *agent)
{
	// Send response back to the agent who made the request
	response = "Yes!  I am the one they call [" + agent->nodeName + ":" + agent->agentName + "].";
	return 0;
}


/// Start the simulation loop
void simulation::start_simulation() {
	//GnuplotPipe gp;
	string node_agent_name = "["+agent->nodeName+":"+agent->agentName+"]";
	cout << node_agent_name << " starting..."<<endl;
	
	// exit with error if unable to start agent
	if(agent->last_error() < 0) {
		cerr<<"error: unable to start "<<node_agent_name
			<<" ("<<agent->last_error()<<") "
			<<cosmos_error_string(agent->last_error())<<endl;
		exit(1);
	} else {
		cout << node_agent_name << " started."<<endl<<endl;
	}

	// Add basic request function
	agent->add_request("are_you_out_there", are_you_out_there, "\n\t\trequest to determine if specific agent exists");

	// Attempt contact with all other agents of the simulation. Exit if not all agents are running.
	if(!all_sim_agents_running()) {
		exit(1);
	}

	// Initialize simulation agents
	init_sim_agents();

	// Simulation loop
	while (agent->running()) {

		cout<<node_agent_name<<" running..."<<endl;

		agent->send_request(agent->find_agent("sat_001", "agent_001", 2.), request, response, 2.);
		if(response.size())	{
			cout<<left<<setw(40)<<"\t[sat_001:agent_001]"<<setw(16)<<"\033[1;32mFOUND\033[0m";
			// ask for their location
			response.clear();
			//agent->send_request(agent->find_agent("sat_001", "agent_001", 2.), "get_position " + time, response, 2.);
			cout<<"\n"<<response<<endl;
		} else {
			cout<<left<<setw(40)<<"\t[sat_001:agent_001]"<<"\033[1;31mNOT FOUND\033[0m"<<endl;
		}
		// Sleep for 5 sec
		COSMOS_SLEEP(5.);
	}

	/*

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
	*/
}


/// Check that all agents in the simulation are running
/**
@return		bool
*/
bool simulation::all_sim_agents_running() {
	cout << "Attempting to find all the agents..." << endl;
	bool out = true;
	response.clear();
	for(int i = 1; i <= num_agents; ++i) {
		string node_name = "sat_" + std::string(3-to_string(i).length(), '0') + to_string(i);
		string agent_name = "agent_" + std::string(3-to_string(i).length(), '0') + to_string(i);
		agent->send_request(agent->find_agent(node_name, agent_name, 2.), request, response, 2.);
		if(response.size())	{
			//cout << "[" << node_name << ":" << agent_name << "] found" << endl;
		} else {
			cout << "Cannot find " << "[" << node_name << ":" << agent_name << "]" << endl;
			out = true;
		}
		response.clear();
	}
	if(out)
		cout << "All agents found!" << endl << endl;
	return out;
}

/// Initialize sim agents
/**
@return n/a
*/
void simulation::init_sim_agents() {
	response.clear();
	// Initialize initial times and states
	std::vector<double> x = {-5,0,5, -5,0,5, -5,0,5};
	std::vector<double> y = { 5,5,5,  0,0,0, -5,-5,-5};
	t0 = std::chrono::duration<double>( std::chrono::system_clock::now().time_since_epoch()).count();
	double t = std::chrono::duration<double>( std::chrono::system_clock::now().time_since_epoch()).count();
	for(int i = 1; i <= num_agents; ++i) {
		string node_name = "sat_" + std::string(3-to_string(i).length(), '0') + to_string(i);
		string agent_name = "agent_" + std::string(3-to_string(i).length(), '0') + to_string(i);
		agent->send_request(agent->find_agent(node_name, agent_name, 2.), request, response, 2.);
		if(response.size())	{
			agent->send_request(agent->find_agent(node_name, agent_name, 2.), "set_initial_time " + to_string(t), response, 2.);
			/*json11::Json json = json11::Json::object {
					{ "x_position", x[i-1] },
					{ "y_position", y[i-1] },
					{ "z_position", 0 },
					{ "x_velocity", 0 },
					{ "y_velocity", 0 },
					{ "z_velocity", 0 }
				};*/
			stringstream ss;
			ss 	<<  "{\"x_position\":" << setprecision(numeric_limits<double>::digits10) << x[i-1] << ","
				<< 	 "\"y_position\":" << setprecision(numeric_limits<double>::digits10) << y[i-1] << ","
				<< 	 "\"z_position\":" << setprecision(numeric_limits<double>::digits10) << 0 << ","
				<< 	 "\"x_velocity\":" << setprecision(numeric_limits<double>::digits10) << 0 << ","
				<< 	 "\"y_velocity\":" << setprecision(numeric_limits<double>::digits10) << 0 << ","
				<< 	 "\"z_velocity\":" << setprecision(numeric_limits<double>::digits10) << 0 << "}";
			agent->send_request(agent->find_agent(node_name, agent_name, 2.), "set_state_vector " + ss.str(), response, 2.);
			this_thread::sleep_for (chrono::milliseconds(10));
			t += 0.01;
		} else {
			std::cerr << "Cannot find " << "[" << node_name << ":" << agent_name << "]" << endl;
			exit(1);
		}
		response.clear();
	}
}

/// Send a request to every agent in the simulation, returns a vector of responses
vector<string> simulation::send_req_to_all_agents(string req) {
	vector<string> all_responses;
	for(int i = 1; i <= num_agents; ++i) {
		string node_name = "sat_" + std::string(3-to_string(i).length(), '0') + to_string(i);
		string agent_name = "agent_" + std::string(3-to_string(i).length(), '0') + to_string(i);
		string response = "";
		agent->send_request(agent->find_agent(node_name, agent_name, 2.), request, response, 2.);
		if(response.size())	{
			agent->send_request(agent->find_agent(node_name, agent_name, 2.), req, response, 2.);
			all_responses.push_back(response);
		} else {
			std::cerr << "Cannot find " << "[" << node_name << ":" << agent_name << "]" << endl;
			exit(1);
		}
	}
	return all_responses;
}













/// Convert a vector into a gnuplot-parsable string
/**
@param	mat			Row/column vector to parse
@param	varname		Variable name to assign the vector. Used by gnuplot.
@return string
*/
string gnuvec(const Eigen::MatrixXd& mat, const string& varname) {
	ostringstream o;
	o << varname << "=\"";
	for(auto it : mat(Eigen::all,Eigen::last))
		o << it << " ";
	o << "\"";
	return o.str();
}

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
void simulation::plot_veh(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y, const vector<double>& trackt, const Eigen::MatrixXd& lx, const Eigen::MatrixXd& obx) {
	gp.sendLine("reset", true);
	gp.sendLine("set title \"Smoothed Particle Hydrodynamics for Agent Control\\n{/*0.85Time = " +to_string(trackt[thead]) + "}\" font \"Arial,16\"", true);
	gp.sendLine("set parametric", true);
	plot_points(x,y);
	plot_lx(lx);
	plot_trails(x,y);
}

/// Plot the SPH particles
/**
@param	x			Matrix of history of x positions of all SPH particles
@param	y			Matrix of history of y positions of all SPH particles
@param	thead		Index position of the head of the vectors
@return n/a
*/
void simulation::plot_points(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y) {
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

/// Plot the loiter circles
/**
@param	lx			Matrix containing [x y] positions of the loiter circles
@return n/a
*/
void simulation::plot_lx(const Eigen::MatrixXd& lx) {
	if(lx.size() != 0) {
		for(auto row : lx.rowwise()) {
			gp.sendLine("set label at " + to_string(row(0)) + "," + to_string(row(1)) + " point pointtype 3 pointsize 2 lt rgb \"goldenrod\"", true);
		}
	}
}

/// Display trail for each particle
/**
@param	x			Matrix of history of x positions of all SPH particles
@param	y			Matrix of history of y positions of all SPH particles
@param	thead		Index position of the head of the vectors
@return n/a
*/
void simulation::plot_trails(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y) {
	for(int i = 0; i < x.rows(); ++i) {
		gp.sendLine(gnutrail(x.row(i),"X"+to_string(i)), true);
		gp.sendLine(gnutrail(y.row(i),"Y"+to_string(i)), true);
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

/// Convert history of particle positions into a gnuplot-parsable line
/**
@param	pos			Matrix of history of positions of SPH particle
@param	varname		Variable name to assign the line. Used by gnuplot.
@param	thead		Index position of the head of the vectors
@return n/a
*/
string simulation::gnutrail(const Eigen::RowVectorXd& pos, const string& varname) {
	ostringstream o;
	o << varname << "=\"";
	for(int i = 0; i < pos.size(); ++i) {
		int tidx = (thead+i+1) % trackMax;
		o << pos(tidx) << " ";
	}
	o << "\"";
	return o.str();
}





/// Get initial time
double simulation::get_initial_time() {
	return t0;
}