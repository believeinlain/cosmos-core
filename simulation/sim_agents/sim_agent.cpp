/********************************************************************
* Copyright (C) 2015 by Interstel Technologies, Inc.
*   and Hawaii Space Flight Laboratory.
*
* This file is part of the COSMOS/core that is the central
* module for COSMOS. For more information on COSMOS go to
* <http://cosmos-project.com>
*
* The COSMOS/core software is licenced under the
* GNU Lesser General Public License (LGPL) version 3 licence.
*
* You should have received a copy of the
* GNU Lesser General Public License
* If not, go to <http://www.gnu.org/licenses/>
*
* COSMOS/core is free software: you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
*
* COSMOS/core is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* Refer to the "licences" folder for further information on the
* condititons and terms to use this software.
********************************************************************/

#include "support/configCosmos.h"
#include "support/elapsedtime.h"
#include "support/json11.hpp"
#include "support/timeutils.h"
#include "agent/agentclass.h"

#include "../sph.h"

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

// Function prototypes
int32_t are_you_out_there(string &request, string &response, Agent *cdata);
int32_t get_initial_time(string &request, string &response, Agent *cdata);
int32_t set_run_state(string &request, string &response, Agent *cdata);
void init_sim_agent();
void HCL(string&);

// ensure the Agent constructor creates only one instance per process
static Agent *agent;
string node_name;
string agent_name;
string node_agent_name;
string request = "are_you_out_there";
string response = "";
bool run = false;
double sleeptime = 5.;
int agent_id;
/// Simulation parameters
param_struct param;
/// Simulation group configurations
group_conf_struct group_conf;
/// Loiter circle locations [x y]
Eigen::MatrixXd lx;
/// Loiter circle radii
Eigen::MatrixXd lR;
/// Reduced density target locations [x y]
Eigen::MatrixXd rdx;
/// Make sure to match group_conf.obs_init
Eigen::MatrixXd obx;

// Simulation parameters

int main(int argc, char **argv)
{
	// Run agent from the command line, with an int id
	if(argc == 2) {
		std::istringstream ss(argv[1]);
		if (!(ss >> agent_id) || !ss.eof()) {
			std::cerr << "Argument " << argv[1] << " is invalid, please provide an int id as argument to this agent." << endl;
			exit(1);
		}
	} else {
		std::cerr << "Please provide an int id as argument to this agent" << endl;
		exit(1);
	}
	node_name = "sat_" + std::string(3-to_string(agent_id).length(), '0') + to_string(agent_id);
	agent_name = "agent_" + std::string(3-to_string(agent_id).length(), '0') + to_string(agent_id);
	node_agent_name = "["+node_name+":"+agent_name+"]";

	// construct agent
	cout << node_agent_name << " starting..."<<endl;
	agent = new Agent(node_name, agent_name, 1.);

	// SPH object
	sph_sim SPH;

	// exit with error if unable to start agent
	if(agent->last_error() < 0) {
		cerr<<"error: unable to start "<<node_agent_name
			<<" ("<<agent->last_error()<<") "
			<<cosmos_error_string(agent->last_error())<<endl;
		exit(1);
	} else {
		cout << node_agent_name << " started."<<endl;
	}

	// add custom request functions for this agent
	agent->add_request("are_you_out_there", are_you_out_there, "\n\t\trequest to determine if specific agent exists");
	agent->add_request("get_initial_time", get_initial_time, "\n\t\trequest to get a simulation parameter, initial time");
	agent->add_request("set_run_state", set_run_state, "\n\t\trequest to set the run state");

	cosmosstruc* c = agent->cinfo;

	init_sim_agent();

	// agent loop
	while (agent->running()) {

		cout<<node_agent_name<<" running..."<<endl;

		// Run if running state is true, set to true by init_sim_agents() in the simulation
		if(run) {
			agent->send_request(agent->find_agent("world", "controller", 2.), request, response, 2.);
			if(response.size())	{
				// Request state vectors
				response.clear();
				agent->send_request(agent->find_agent("world", "controller", 2.), "get_state_vectors", response, 2.);

				// Use pseudo-HCL to time-align state
				HCL(response);

				// Update sph state vector to current positions
				//SPH.sph_update_state(agent->cinfo->get_json<vector<statestruct>>("state"), agent_id);

				// Calculate next waypoint via SPH
				//SPH.sph_sim_step(rdx,lx,lR);

				// Send world controller updated agent state
				//agent->send_request(agent->find_agent("world", "controller", 2.), "send_world_new_state " + agent->cinfo->get_json<statestruct>("state["+to_string(agent_id-1)+"]"), response, 2.);


				cout<<left<<setw(40)<<"\t[world:controller]"<<setw(16)<<"\033[1;32mFOUND\033[0m";
				// ask for their location
				response.clear();
			} else {
				cout<<left<<setw(40)<<"\t[world:controller]"<<"\033[1;31mNOT FOUND\033[0m"<<endl;
				run = false;
				sleeptime = 5.;
			}
		}
		
		// Sleep for 5 sec
		COSMOS_SLEEP(sleeptime);
	}
	return 0;
}


/// Basic role call function
int32_t are_you_out_there(string & request, string &response, Agent *) {
	// Send response back to the agent who made the request
	response = "Yes!  I am the one they call " + node_agent_name + ".";
	return 0;
}

/// Get simulation initial time
int32_t get_initial_time(string &request, string &response, Agent *agent) {
	cout<<"\tincoming request          = <"<<request<<">"<<endl;
	response = to_string(agent->cinfo->get_value<double>("state["+to_string(agent_id-1)+"].timestamp"));
	return 0;
}

/// Set run state
int32_t set_run_state(string &request, string &response, Agent *) {
	cout<<"\tincoming request          = <"<<request<<">"<<endl;

	// remove function call and space
	request.erase(0,14);

	// read in request arg
	istringstream ss(request);
	ss>>std::boolalpha>>run;
	// how often to run the sph update, for our purposes doesn't necessarily have to be equal to dt
	sleeptime = 0.5;
	return 0;
}

void init_sim_agent() {
	//param = SPH.get_param();
	//group_conf = SPH.get_group_conf();
	//group_conf.num_veh(0) = 9;
	// Loiter circle position
	lx.resize(1,2);
	lx << 28,0;
	lR.resize(1,1);
	lR << 5;
	// Attractors
	rdx.resize(1,2);
	rdx << 28,0;
	// Obstacle positions
	// Make sure to match group_conf.obs_init
	obx.resize(5,2);
	obx << 	 7, 0,
			12, 4,
			16, 2,
			 9,-2,
			22,-6;
}

void HCL(string &state) {
	// parse json state
	string error;
	json11::Json parsed = json11::Json::parse(state,error);
	if(error.empty()) {
		// multiply velocity uvw by time difference between the current time and timestamp, then add to position xyz
		double t = std::chrono::duration<double>( std::chrono::system_clock::now().time_since_epoch()).count();
		for(size_t i = 0; i < parsed["state"].array_items().size(); ++i) {
			if(int(i+1) != agent_id) {
				double dt = parsed["state"][i]["timestamp"].number_value() - t;
				agent->cinfo->set_value<double>("state["+to_string(i)+"].x_position", parsed["state"][i]["x_position"].number_value() + parsed["state"][i]["x_velocity"].number_value() * dt);
				agent->cinfo->set_value<double>("state["+to_string(i)+"].y_position", parsed["state"][i]["y_position"].number_value() + parsed["state"][i]["y_velocity"].number_value() * dt);
				agent->cinfo->set_value<double>("state["+to_string(i)+"].z_position", parsed["state"][i]["z_position"].number_value() + parsed["state"][i]["z_velocity"].number_value() * dt);
			}
		}
		cout << agent->cinfo->get_json<statestruct>("state") << endl;
	}
}



// Notes to self:
// Important to declare SPH object after new Agent is called, otherwise segfaults. Probably has to do with
// how agents need to have their own process id or something. Needs to be the first object created in the executable, or something.