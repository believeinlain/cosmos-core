#ifndef SPH_STRUCTS_H
#define SPH_STRUCTS_H

#include "Eigen/Dense"

//using namespace Eigen;
using namespace std;

enum particle_type_enum {veh, obs, rd };

/// A 3D matrix, N x M x 3
struct Matrix3D {
	Eigen::MatrixXd z0;
	Eigen::MatrixXd z1;
	Eigen::MatrixXd z2;
	Matrix3D()  {}
	Matrix3D(Eigen::MatrixXd a, Eigen::MatrixXd b, Eigen::MatrixXd c) : z0(a), z1(b), z2(c) {}
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
	Eigen::VectorXd x, y, z;

	// initial velocities
	Eigen::VectorXd u, v, w;
};

/// Limits for speed and acceleration
struct veh_limits_struct {
	Eigen::VectorXd vmin, vmax;
	Eigen::VectorXd turning_radius;
};

/// Positions
struct pos_struct {
	Eigen::VectorXd x, y, z;
};

/// Structure containing all the properties of each SPH particle
struct prop_struct {
	/// Minimum velocity constraint
	Eigen::MatrixXd vmin;
	/// Maximum velocity constraint
	Eigen::MatrixXd vmax;
	/// Turning radius constraint
	Eigen::MatrixXd turning_radius;
	/// Maximum acceleration constraint
	Eigen::MatrixXd amax;
	/// Kernel width
	Eigen::MatrixXd h;
	/// Mass
	Eigen::MatrixXd m;
	/// Viscosity
	Eigen::MatrixXd mu;
	/// Bulk modulus
	Eigen::MatrixXd K;
	/// Group number
	Eigen::MatrixXd group;
	/// Particle type (veh, obs, or rd)
	Eigen::MatrixXd particle_type;
	/// h_ij matrix
	Eigen::MatrixXd hij;
	/// Kernel type (1 or 2)
	Eigen::MatrixXd kernel_type;
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
	Eigen::VectorXi num_veh;

	/// Struct of vectors of initial positions and velocities for the vehicle groups
	pos_vel_struct veh_init;

	/// Vector of smoothing widths for each group
	Eigen::VectorXd veh_h;

	/// Limits for speed and acceleration
	veh_limits_struct veh_limits;

	/// Total number of obstacle particles
	int num_obs;

	/// Vector of size of particles (times 0.5)
	Eigen::VectorXd obs_h;

	/// Vector of positions of the obstacles
	pos_struct obs_init;

	/// Total number of reduced density particles
	int num_rd;

	/// Which group does each reduced density particle belong to?\n
	/// Group number corresponds to array index for num_veh\n
	/// -1 means not active
	Eigen::VectorXd rd_group;

	/// Vector of initial positions and velocities for the reduced particles
	pos_vel_struct rd_init;

	/// Smoothing width for each reduced particle group
	Eigen::VectorXd rd_h;

	/// Total number of loiter circles
	int num_loiter;

	/// Which group does each loiter circle belong to?\n
	/// Group number corresponds to array index for num_veh\n
	/// -1 means not active
	Eigen::VectorXi loiter_group;

};


#endif