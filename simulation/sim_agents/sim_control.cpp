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
#include <string>

// The request function prototype
int32_t are_you_out_there(string &request, string &response, Agent *cdata);

// ensure the Agent constructor creates only one instance per process
static Agent *agent;
string node_name = "world";
string agent_name = "control";
string node_agent_name = "["+node_name+":"+agent_name+"]";

string request = "are_you_out_there";
string response = "";

int main(int argc, char **argv)
{
	/*
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

	cosmosstruc* c = agent->cinfo;

	// agent loop
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
		agent->send_request(agent->find_agent("sat_002", "agent_002", 2.), request, response, 2.);
		if(response.size())	{
			cout<<left<<setw(40)<<"\t[sat_002:agent_002]"<<setw(16)<<"\033[1;32mFOUND\033[0m";
			// ask for their location
			response.clear();
			//agent->send_request(agent->find_agent("sat_002", "agent_002", 2.), "get_position " + time, response, 2.);
			cout<<"\n"<<response<<endl;
		} else {
			cout<<left<<setw(40)<<"\t[sat_002:agent_002]"<<"\033[1;31mNOT FOUND\033[0m"<<endl;
		}
		
		// Sleep for 5 sec
		COSMOS_SLEEP(5.);
	}*/
	return 0;
}


bool agents_are_running() {

	return false;
}


int32_t are_you_out_there(string & request, string &response, Agent *)
{
	// Send response back to the agent who made the request
	response = "Yes!  I am the one they call " + node_agent_name + ".";
	return 0;
}