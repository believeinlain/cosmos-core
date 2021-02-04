#include <iostream>
 
#include "simulation.h"
#include "agent/agentclass.h"

using namespace std;
static Agent* agent;


int main()
{
	// Create simulation control agent
	string node_name = "world"; 
	string agent_name = "controller";
	string node_agent_name = "["+node_name+":"+agent_name+"]";
	agent = new Agent(node_name, agent_name, 1.);
	
	// Simulation class
	simulation sim;
	cout << "Starting SPH simulation..." << endl;
	sim.start_simulation(agent);

	cout << "Ending SPH simulation." << endl;
}