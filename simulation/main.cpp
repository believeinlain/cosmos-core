#include <iostream>
 
#include "simulation.h"

using namespace std;

int main()
{
	simulation sim;
	cout << "Starting SPH simulation..." << endl;
	sim.start_simulation();

	cout << "Ending SPH simulation." << endl;
}
