#ifndef SPH_SIM_H
#define SPH_SIM_H

#include <math.h>

#include "matrix.h"
#include "sph_structs.h"


class sph_sim {
private:
	// Loiter circle x,y,R
	double lx, lR;

	// Time
	double t, t0;

	// List of properties for each SPH particle
	prop_struct prop;

	// Parameters for the SPH equations
	param_struct param;

	// Group configuration (number of vehicles, vehicle constraints, etc)
	group_conf_struct group_conf;

	// Particle states
	// A matrix containing the [x y z] positions and [u v w] of all the SPH particles. 
	// Each particle is stored in one row so the format is:
	// 			[ x1 y1 z1 u1 v1 w1 ]
	// states = [ x2 y2 z2 u2 v2 w2 ]
	// 			[        ...        ]
	matrix states;

	// Total number of SPH particles
	unsigned int npart;

	// Number of vehicles
	unsigned int nveh, nobs, nrd;

	// Reference density
	double rho0;

	// Initialize properties for each SPH particle
	void init_prop();

	// Compute h_ij matrix
	void compute_hij();

	// Create a matrix kernel_type that tells which kernel to use.
	// 1 is for vehicle-reduced density particle interactions,
	// 2 is for all others
	void kernel_type();

	// Initialize positions/velocities
	void init_states();

	// 2D initialization, use a 2D hexagonal (close packed) lattice
	void init2d();

	// 3D initialization, use a 3D hexagonal (close packed) lattice
	void init3d();
	
public:
	// Class constructor
	// param		structure containing the SPH parameters
	// group_conf	structure containing the SPH group configuration
	// t0			option, if t0 is not given it is set to 0
	sph_sim(param_struct param, group_conf_struct group_conf, double t0 = 0);

};

#endif