#ifndef SPH_STRUCTS_H
#define SPH_STRUCTS_H

#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

enum particle_type_enum {veh, obs, rd };

/// A 3D matrix, N x M x 3
struct Matrix3D {
	MatrixXd z0;
	MatrixXd z1;
	MatrixXd z2;
	Matrix3D()  {}
	Matrix3D(MatrixXd a, MatrixXd b, MatrixXd c) : z0(a), z1(b), z2(c) {}
};

/// Gain coefficients
struct gain_struct {
	/// Gain coefficient for the SPH forces
	double sph;
	/// Gain coefficient for the external force
	double ext;
	/// Gain coefficient for the drag force
	double drag;
};

/// Scales accel due to vehicles/obstacles/reduced density particles
struct accel_struct {
	/// Scaling constant for SPH vehicle forces
	double veh;
	/// Scaling constant for SPH obstacle forces
	double obs;
	/// Scaling constant for SPH attractor forces
	double rd;
};

/// Initial positions and velocities for the vehicle groupps
struct pos_vel_struct {
	// initial positions
	VectorXd x, y, z;

	// initial velocities
	VectorXd u, v, w;
};

/// Limits for speed and acceleration
struct veh_limits_struct {
	VectorXd vmin, vmax;
	VectorXd turning_radius;
};

/// Positions
struct pos_struct {
	VectorXd x, y, z;
};

/// Structure containing all the properties of each SPH particle
struct prop_struct {
	/// Minimum velocity constraint
	MatrixXd vmin;
	/// Maximum velocity constraint
	MatrixXd vmax;
	/// Turning radius constraint
	MatrixXd turning_radius;
	/// Maximum acceleration constraint
	MatrixXd amax;
	/// Kernel width
	MatrixXd h;
	/// Mass
	MatrixXd m;
	/// Viscosity
	MatrixXd mu;
	/// Bulk modulus
	MatrixXd K;
	/// Group number
	MatrixXd group;
	/// Particle type (veh, obs, or rd)
	MatrixXd particle_type;
	/// h_ij matrix
	MatrixXd hij;
	/// Kernel type (1 or 2)
	MatrixXd kernel_type;
};


/// SPH simulation parameters
struct param_struct {
	/// Dimension of the simulation (2 or 3)
	int ndim = 2;

	/// Gain coefficients
	gain_struct gain;

	/// Scaling constant for particles
	accel_struct accel;

	/// Reynolds number
	double Re;

	/// Timestep for the SPH simulation
	double dt;
};

/// The groups of vehicle, obstacles, and reduced density particles
struct group_conf_struct {
	/// Vector containing the number of vehicles in each group
	VectorXi num_veh;

	/// Struct of vectors of initial positions and velocities for the vehicle groups
	pos_vel_struct veh_init;

	/// Vector of smoothing widths for each group
	VectorXd veh_h;

	/// Limits for speed and acceleration
	veh_limits_struct veh_limits;

	/// Total number of obstacle particles
	int num_obs;

	/// Vector of size of particles (times 0.5)
	VectorXd obs_h;

	/// Vector of positions of the obstacles
	pos_struct obs_init;

	/// Total number of reduced density particles
	int num_rd;

	/// Which group does each reduced density particle belong to?\n
	/// Group number corresponds to array index for num_veh\n
	/// -1 means not active
	VectorXd rd_group;

	/// Vector of initial positions and velocities for the reduced particles
	pos_vel_struct rd_init;

	/// Smoothing width for each reduced particle group
	VectorXd rd_h;

	/// Total number of loiter circles
	int num_loiter;

	/// Which group does each loiter circle belong to?\n
	/// Group number corresponds to array index for num_veh\n
	/// -1 means not active
	VectorXi loiter_group;

};


#endif