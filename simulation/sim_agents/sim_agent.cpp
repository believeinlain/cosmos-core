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
#include "support/timeutils.h"
#include "agent/agentclass.h"

#include <iostream>
#include <sstream>
#include <string>

// Function prototypes
int32_t are_you_out_there(string &request, string &response, Agent *cdata);
int32_t get_initial_time(string &request, string &response, Agent *cdata);
int32_t set_initial_time(string &request, string &response, Agent *cdata);

// ensure the Agent constructor creates only one instance per process
static Agent *agent;
string node_name;
string agent_name;
string node_agent_name;
string request = "are_you_out_there";
string response = "";

// Simulation parameters
/// Initial time
double t0;

int main(int argc, char **argv)
{
	int agent_id;

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
	agent->add_request("set_initial_time", set_initial_time, "\n\t\trequest to set a simulation parameter, initial time");

	cosmosstruc* c = agent->cinfo;

	// agent loop
	while (agent->running()) {

		cout<<node_agent_name<<" running..."<<endl;

		agent->send_request(agent->find_agent("world", "controller", 2.), request, response, 2.);
		if(response.size())	{
			cout<<left<<setw(40)<<"\t[world:controller]"<<setw(16)<<"\033[1;32mFOUND\033[0m";
			// ask for their location
			response.clear();
			//agent->send_request(agent->find_agent("sat_002", "agent_002", 2.), "get_position " + time, response, 2.);
			cout<<"\n"<<response<<endl;
		} else {
			cout<<left<<setw(40)<<"\t[world:controller]"<<"\033[1;31mNOT FOUND\033[0m"<<endl;
		}
		// Sleep for 5 sec
		COSMOS_SLEEP(5.);
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
	response = to_string(t0);
	return 0;
}

/// Set simulation initial time
int32_t set_initial_time(string & request, string &response, Agent *) {
	cout<<"\tincoming request          = <"<<request<<">"<<endl;

	// remove function call and space
	request.erase(0,17);
	
	// read in mjdtime
	stringstream ss;
	ss<<request;
	ss>>t0;
	cout << "Initial time t0 set to: " << to_string(t0) << endl;


	return 0;
}