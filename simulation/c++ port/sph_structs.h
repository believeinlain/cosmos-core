#ifndef SPH_STRUCTS_H
#define SPH_STRUCTS_H

#include <vector>

#include "matrix.h"

using namespace std;

enum particle_type_enum {veh, obs, rd };

// Gain coefficients
struct gain_struct {
	double sph;
	double ext;
	double drag;
};

// Scales accel due to vehicles/obstacles/reduced density particles
struct accel_struct {
	double veh;
	double obs;
	double rd;
};

// Initial positions and velocities for the vehicle groupps
struct pos_vel_struct {
	// initial positions
	vector<double> x, y, z;

	// initial velocities
	vector<double> u, v, w;
};

// Limits for speed and acceleration
struct veh_limits_struct {
	matrix vmin, vmax;
	matrix turning_radius;
};

// Positions
struct pos_struct {
	vector<double> x, y, z;
};

// SPH particle properties
// Structure containing all the properties of each SPH
struct prop_struct {
	matrix vmin;					// Minimum velocity constraint
	matrix vmax;					// Maximum velocity constraint
	matrix turning_radius;			// Turning radius constraint
	matrix amax;					// Maximum acceleration constraint
	matrix h;						// Kernel width
	matrix m;						// Mass
	matrix mu;						// Viscosity
	matrix K;						// Bulk modulus
	matrix group;					// Group number
	matrix particle_type;		// Particle type (veh, obs, or rd)
	matrix hij;						// h_ij matrix
	matrix kernel_type;
};


// SPH simulation parameters
struct param_struct {
	// Dimension of the simulation (2 or 3)
	unsigned int ndim = 2;

	gain_struct gain;
	accel_struct accel;

	// Reynolds number
	double Re;

	// Timestep for the SPH simulation
	double dt;
};

// The groups of vehicle, obstacles, and reduced density particles
struct group_conf_struct {
	// Vector containing the number of vehicles in each group
	vector<unsigned int> num_veh;

	// Struct of vectors of initial positions and velocities for the vehicle groups
	pos_vel_struct veh_init;

	// Vector of smoothing widths for each group
	matrix veh_h;

	// Limits for speed and acceleration
	veh_limits_struct veh_limits;

	// Total number of obstacle particles
	unsigned int num_obs;

	// Vector of size of particles (times 0.5)
	matrix obs_h;

	// Vector of positions of the obstacles
	pos_struct obs_init;

	// Total number of reduced density particles
	unsigned int num_rd;

	// Which group does each reduced density particle belong to?
	// Group number corresponds to array index for num_veh
	// -1 means not active
	matrix rd_group;

	// Vector of initial positions and velocities for the reduced particles
	pos_vel_struct rd_init;

	// Smoothing width for each reduced particle group
	matrix rd_h;

	// Total number of loiter circles
	unsigned int num_loiter;

	// Which group does each loiter circle belong to?
	// Group number corresponds to array index for num_veh
	// -1 means not active
	// loiter_group

};


#endif