#ifndef SPH_SIM_H
#define SPH_SIM_H

#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "sph_structs.h"

using namespace Eigen;
using namespace std;


double kernel(double r, double h, int type);
double kernel_grad(double r, double h, int type);
MatrixXi find(const MatrixXd& A);
MatrixXd index(const MatrixXd& m, const MatrixXi& I);
MatrixXd vseq(int val0, int val_last);
MatrixXi sort(const MatrixXd& c);
MatrixXd append_right(const MatrixXd& m, const MatrixXd& app);
MatrixXd append_down(const MatrixXd& m, const MatrixXd& app);

class sph_sim {
private:
	// Loiter circle x,y,z, R
	MatrixXd lx;
	MatrixXd lR;

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
	MatrixXd states;

	// Total number of SPH particles
	int npart;

	// Number of vehicles
	int nveh, nobs, nrd;

	// Reference density
	const double rho0 = 1;

	// Default simulation parameters
	void init();

	// Initialize properties for each SPH particle
	void init_prop();

	// Resize vectors of prop, adds n more rows. For internal use
	void resize_prop(int rows);

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

	// Return the right hand side of the SPH momentum equation
	MatrixXd sph_rhs();

	// Compute the distance, vector, and unit vector between particles i and j
	tuple<MatrixXd, Matrix3D, Matrix3D> sph_compute_dij();

	// Compute the masking function and the indices that are non-zero
	//			{ 0 if dij>2*hij
	// M(i,j) = { 0 if particle j is not a vehicle and group(i)~=group(j)
	//			{ 0 if particle i is not a vehicle (except M(i,i)=1)
	//			{ 1 else
	tuple<MatrixXd, MatrixXi> sph_compute_mask(const MatrixXd& dij);

	MatrixXd sph_compute_density(const MatrixXd& dij, const MatrixXd& Mask, const MatrixXi& MaskI);

	// Equation of state to compute the pressure
	MatrixXd sph_compute_pressure(const MatrixXd& rho);

	// Compute the viscous forces
	MatrixXd sph_compute_pi(const MatrixXd& rho, const MatrixXd& dij, const Matrix3D& rij, const Matrix3D& unit_ij,
								const MatrixXd& gradW, const MatrixXd& Mask, const MatrixXi& MaskI);
	
	// Compute the external force on vehicles to drive them toward a loiter circle
	tuple<MatrixXd,MatrixXd,MatrixXd> external_force();

	// Compute the rate of change of SPH.states, i.e., the velocity
	// and accelerations, while applying vehicle constraints
	MatrixXd sph_compute_rates(const MatrixXd& DvDt);

	// Apply velocity constraints
	void constrain_vel();

	
public:
	// Class constructor
	// param		structure containing the SPH parameters
	// group_conf	structure containing the SPH group configuration
	// t0			option, if t0 is not given it is set to 0
	sph_sim();
	sph_sim(param_struct param, group_conf_struct group_conf, double t0 = 0);

	// Update or change the SPH properties to match the properties in the arguments param and group_conf
	void sph_update_properties(const param_struct& param, const group_conf_struct& group);

	// Take a single time-step forward in the simulation
	void sph_sim_step(MatrixXd rdx, MatrixXd lx, MatrixXd lR);

	// GETTERS
	// Return the current time in the SPH simulation
	double get_time();
	// Return the initial time for the SPH simulation
	double get_initial_time();
	// Return the time step to be used for the SPH simulation
	double get_dt();

	// Return a matrix containing the [x y z] positions and [u v w] velocities of all SPH particles.
	// Each particle is stored in one row:
	//			[ x0 y0 z0 u0 v0 w0 ]
	// states = [ x1 y1 z1 u1 v1 w1 ]
	//			[		 ...		]
	MatrixXd get_states();

	// Return the total number of particles in the simulation
	int get_npart();
	// Return the number of vehicles in the simulation
	int get_nveh();
	// Return the number of obstacles in the simulation
	int get_nobs();
	// Return the number of reduced density (attractor) particles in the simulation
	int get_nrd();

	// Return a column vector containing x/y/z positions or u/v/w velocities of all the SPH particles
	MatrixXd get_x();
	MatrixXd get_y();
	MatrixXd get_z();
	MatrixXd get_u();
	MatrixXd get_v();
	MatrixXd get_w();

	// Returns a prop_struct containing all the properties of each SPH particle in the simulation
	// Members of prop_struct:
	// MatrixXd vmin				Minimum velocity constraint
	// MatrixXd vmax				Maximum velocity constraint
	// MatrixXd turning_radius		Turning radius constraint
	// MatrixXd amax				Maximum acceleration constraint
	// MatrixXd h					Kernel width
	// MatrixXd m					Mass
	// MatrixXd mu					Viscosity
	// MatrixXd K					Bulk modulus
	// MatrixXd group				Group number
	// MatrixXd particle_type		Particle type (veh, obs, or rd)
	// MatrixXd hij					h_ij matrix
	// MatrixXd kernel_type			kernel type
	prop_struct get_prop();

	// Returns a param_struct containing all the parameters used in the SPH simulation
	// Members of param_struct:
	// int param.ndim				dimension of the simulation (2 or 3)
	// double param.gain.sph		gain coefficient for the SPH forces
	// double param.gain.ext		gain coefficient for the external force
	// double param.gain.drag		gain coefficient for the drag force
	// double param.accel.veh		scaling constant for SPH vehicle forces
	// double param.accel.obs		scaling constant for SPH obstacle forces
	// double param.accel.rd		scaling constant for SPH attractor forces
	// double param.Re				Reynolds number
	// double param.dt				Time step
	param_struct get_param();

	// Return a group_conf_struct containing the group configuration
	group_conf_struct get_group_conf();

};

#endif