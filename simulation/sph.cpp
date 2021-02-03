#include "sph.h"

/// Evaluates the smoothing kernel function for the SPH equations
/**
The equations:

<center><table border="0"><tr>
<td style="text-align:left; padding-right:15mm">
\f$
\begin{align*}
W_1(\vec{r_{i j}},h_{i j})=\frac{C_1}{h_{i j}^d}\left\{\begin{matrix}
	4-6s^2+3s^3	& if & 0 \leq  s \leq 1 \\
	(2-s)^3		& if & 1 \leq s \leq 2 \\
	0			& if & \hfill s > 2
\end{matrix}\right.
\end{align*}
\f$
</td><td>
and
</td><td style="text-align:right; padding-left:15mm">
\f$
\begin{align*}
W_2(\vec{r_{i j}},h_{i j})=\frac{C_2}{h_{i j}^d}\left\{\begin{matrix}
	s^2-4s+4	& if & s \leq 2 \\
	0			& if & s > 2
\end{matrix}\right.
\end{align*}
\f$
</td>
</tr></table></center>

are the two kernel functions used in this simulation, for type == 1 and type == 2, respectively.

@param	r		relative distance between particles i and j
@param	h		smoothing width, \f$ max(h_i,h_j) \f$ if they differ.
@param	type	int type value for picking the kernel function to use. Type 1 for attracting forces and type 2 for repelling vehicle-vehicle interactions.

@return W		double value result of the kernel function
*/
double kernel(double r, double h, int type) {
	double s = r/h;
	
	double W = 0;

	// Cubic spline kernel, used for vehicle-reduced density (type 1) particle interactions.
	W = W + (type==1)*( ( 1.0-3.0/2.0*pow(s,2.0) + 3.0/4.0*pow(s,3.0) ) * (s<=1) + ( 1.0/4.0*pow((2.0-s),3.0) ) * ( (s > 1.0) * (s <= 2.0) ) ) / (M_PI*pow(h,3.0));

	// Second-order spline kernel for vehicle-vehicle interactions. Generates repelling forces if two particles are too close. 
	W = W + (type==2)*15.0/(16.0*M_PI*pow(h,3.0))*(pow(s,2.0)/4.0-s+1.0)*(s <= 2.0);
	return W;
}


/// Evaluates the derivative of the smoothing kernel function
/**
The gradients of the two kernel functions are:

<center><table border="0"><tr>
<td style="text-align:left; padding-right:15mm">
\f$
\begin{align*}
\frac{d W_1}{d r}(\vec{r_{i j}},h_{i j})=\frac{C_1}{h_{i j}^{d+1}}\left\{\begin{matrix}
	-12s+9s^2	& if & 0 \leq  s \leq 1 \\
	-3(2-s)^2	& if & 1 \leq s \leq 2 \\
	0			& if & \hfill s > 2
\end{matrix}\right.
\end{align*}
\f$
</td><td>
and
</td><td style="text-align:right; padding-left:15mm">
\f$
\begin{align*}
\frac{d W_2}{d r}(\vec{r_{i j}},h_{i j})=\frac{C_2}{h_{i j}^{d+1}}\left\{\begin{matrix}
	2s-4	& if & s \leq 2 \\
	0		& if & s > 2
\end{matrix}\right.
\end{align*}
\f$
</td>
</tr></table></center>

for type == 1 and type == 2, respectively.

Note that this returns a scalar value which is the magnitude of the gradient W.

The direction (if needed) must be computed separately.

@param	r		relative distance between particles i and j
@param	h		smoothing width, \f$ max(h_i,h_j) \f$ if they differ.
@param	type	int type value for picking the kernel function to use. Type 1 for attracting forces and type 2 for repelling vehicle-vehicle interactions.

@return W		double value result of the kernel function
*/
double kernel_grad(double r, double h, int type) {
	double s = r/h;

	double dWdr = 0;

	// Cubic spline kernel, used for vehicle-reduced density (type 1) particle interactions.
	dWdr = dWdr + (type==1)*( ( -3.0*s + 9.0/4.0*pow(s,2.0))*(s<1) + ( -3.0/4.0*pow((2.0-s),2.0) )*( (s>=1.0)*(s<=2.0)  ) )/(M_PI*pow(h,4.0));

	// Quadratic kernel, used for all other (type 2) interactions
	dWdr = dWdr + (type==2)*15.0/(16.0*M_PI*pow(h,4.0))*(s/2.0-1.0)*(s<=2.0);

	return dWdr;
}

/// Find indices of nonzero elements of a matrix
/**
@param	m		Matrix to find indices of.

@return idxs	Vector of indices of nonzero elements. As row if A is a row vector, as a column vector otherwise.
*/
Eigen::MatrixXi find(const Eigen::MatrixXd& m) {
	Eigen::MatrixXi idxs(1,m.size());
	int size = 0;
	for(Eigen::Index i=0; i<m.size(); ++i) {
		if(m(i)) {
			idxs(0,size) = i;
			size++;
		}
	}
	if(size > 0) {
		idxs.conservativeResize(1,size);
		if(m.rows() > 1) {
			idxs.transposeInPlace();
		}
	} else {
		idxs.resize(0,0);
	}
	
	return idxs;
}

/// Generate a row vector with sequentially increasing values
/**
Returns a row vector of doubles in rising sequence by 1 (eg: {0, 1, 2, 3, ...}), inclusive.

@param	val0	int value to start at
@param	valn	int value to end at. Inclusive.

@return rseq	Row vector
*/
Eigen::MatrixXd vseq(int val0, int valn) {
	Eigen::RowVectorXd rseq;
	rseq.setLinSpaced(valn-val0+1, val0, valn);
	return rseq.matrix();
}

/// Get sorted indices of a matrix in ascending order
/**
@param	m		Matrix to find sorted indices for

@return sorted	Row vector of sorted indices
*/
Eigen::MatrixXi sort(const Eigen::MatrixXd& m) {
	Eigen::MatrixXi sorted(1,m.size());
	for(int i = 0; i < m.size(); ++i) {
		sorted(i) = i;
	}
	sort(sorted.data(), sorted.data()+sorted.size(),[&](int i, int j) {
		return m(i) < m(j);
	});
	return sorted;
}

/// Append matrices together horizontally
/**
Append a matrix app to the right of matrix m. Note that the number of rows must match. Argument matrices are not modified.

@param	m		Matrix to find append onto
@param	app		Matrix to append

@return out		Appended matrix
*/
Eigen::MatrixXd append_right(const Eigen::MatrixXd& m, const Eigen::MatrixXd& app) {
	int off = m.cols();
	Eigen::MatrixXd out(m.rows(), off + app.cols());
	out.topLeftCorner(m.rows(), m.cols()) = m;
	out(Eigen::seq(0,Eigen::last),Eigen::seq(off,Eigen::last)) = app;
	return out;
}

/// Append matrices together vertically
/**
Append a matrix app below matrix m. Note that the number of columns must match. Argument matrices are not modified.

@param	m		Matrix to find append onto
@param	app		Matrix to append

@return out		Appended matrix
*/
Eigen::MatrixXd append_down(const Eigen::MatrixXd& m, const Eigen::MatrixXd& app) {
	int off = m.rows();
	Eigen::MatrixXd out(off + app.rows(), m.cols());
	out.topLeftCorner(m.rows(), m.cols()) = m;
	out(Eigen::seq(off,Eigen::last),Eigen::seq(0,Eigen::last)) = app;
	return out;
}



/// Constructor
/**

*/
sph_sim::sph_sim() {
	// Initialize default settings
	init();

	// Setup SPH properties
	init_prop();
	compute_hij();
	kernel_type();

	// Initialize positions and velocities
	init_states();
}

/// Update SPH properties
/**
Update or change the SPH properties to match the properties in the arguments param and group_conf. For use in the simulation.

@param	param	Updated param_struct
@param	group	Updated group_conf_struct

@return n/a
*/
void sph_sim::sph_update_properties(const param_struct& param, const group_conf_struct& group) {
	// Define parameters
	this->param = param;

	// Define group configuration
	this->group_conf = group_conf;

	// Setup SPH properties
	init_prop();
	compute_hij();
}

/// Take a single time-step forward in the simulation
/**
Advance the simulation time clock. Fourth-order Runge-Kutta method is used to approximate optimal path (?).

@param	rdx		Updated reduced density particle positions
@param	lx		Updated loiter circle positions
@param	lR		Updated loiter circle radii

@return n/a
*/
void sph_sim::sph_sim_step(const Eigen::MatrixXd& rdx,const Eigen::MatrixXd& lx,const Eigen::MatrixXd& lR) {
	// NOTE: this simulation only uses the first 3 args, 2 more args are in the octave code
	// Set reduced density particle positions
	if(rdx.size() != 0 && group_conf.num_rd > 0) {
		if(rdx.rows() != group_conf.num_rd && rdx.rows() != 1) {
			cout << "Error in input to sph_sim_step:" << endl;
			cout << "Number of reduced density particle positions must equal the number of reduced density particles." << endl;
			cout << "rdx.rows(): " << rdx.rows() << "; group_conf.num_rd: " << group_conf.num_rd << ";" << endl;
			throw "Error in input to sph_sim_step";
		}
		int I1 = nveh + group_conf.num_obs;
		states.col(0)(Eigen::seq(I1,Eigen::last)) = rdx.col(0);
		states.col(1)(Eigen::seq(I1,Eigen::last)) = rdx.col(1);
		if(param.ndim == 3) {
			states.col(2)(Eigen::seq(I1,Eigen::last)) = rdx.col(2);
		}
	}

	// Set the loiter circle positions and radii
	if(lx.size() != 0 && group_conf.num_loiter > 0) {
		if(lx.rows() != group_conf.num_loiter) {
			cout << "Error in input to sph_sim_step:" << endl;
			cout << "Number of loiter circle positions must equal the number of loiter circles." << endl;
			cout << "lx.rows(): " << lx.rows() << "; group_conf.num_loiter: " << group_conf.num_loiter << ";" << endl;
			throw "Error in input to sph_sim_step";
		}
		this->lx = lx;

		// Set radius
		if(lR.size() != 0) {
			this->lR = lR;
		} else {
			// Radius = minimum turning radius
			this->lR.resize(group_conf.num_loiter,1);
			for(int i = 0; i < group_conf.num_loiter; ++i) {
				int group_num = group_conf.loiter_group(i);
				this->lR(i) = max( 	0.5 * group_conf.veh_h(group_num) / sin( M_PI / max( 2 , group_conf.num_veh(group_num) ) ),
									group_conf.veh_limits.turning_radius(group_num) );
			}
		}
	}

	// Assume the background velocity is changing slower than the SPH timestep so just load it once
	// if exist(velfunc) //NOTE: velfunc arg not used
	Eigen::MatrixXd uvw = Eigen::MatrixXd::Zero(states.rows(),3);

	// Fourth order Runge-Kutta timestep
	Eigen::MatrixXd k1 = sph_rhs();
	k1(Eigen::all,Eigen::seq(0,2)) = k1(Eigen::all,Eigen::seq(0,2)) + uvw;
	
	sph_sim tmp = *this;
	tmp.states = tmp.states + tmp.param.dt/2.0*k1;
	Eigen::MatrixXd k2 = tmp.sph_rhs();
	k2(Eigen::all,Eigen::seq(0,2)) = k2(Eigen::all,Eigen::seq(0,2)) + uvw;

	tmp = *this;
	tmp.states = tmp.states + tmp.param.dt/2.0*k2;
	Eigen::MatrixXd k3 = tmp.sph_rhs();
	k3(Eigen::all,Eigen::seq(0,2)) = k3(Eigen::all,Eigen::seq(0,2)) + uvw;

	tmp = *this;
	tmp.states = tmp.states + tmp.param.dt*k3;
	Eigen::MatrixXd k4 = tmp.sph_rhs();
	k4(Eigen::all,Eigen::seq(0,2)) = k4(Eigen::all,Eigen::seq(0,2)) + uvw;

	states = states + param.dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
	// Increment time
	t = t + param.dt;

	// Constrain the velocity
	constrain_vel();
}

/// Initialize values for demo simulation
/**
@return n/a
*/
void sph_sim::init() {
	// Initialize SPH simulation parameters
	rho0 = 1;
	t0 = 0;
	t = t0;
	param.ndim = 2;
	param.gain.sph = 1.0;
	param.gain.ext = 0.25;
	param.gain.drag = 10.0;
	param.accel.veh = 1.0;
	param.accel.obs = 1.0;
	param.accel.rd = 0.25;
	param.Re = 10.0;
	param.dt = 0.01;

	// Initialize group configurations
	// The groups of vehicles, obstacles, and reduced density particles
	// (Dimensional parameters)
	
	// Number of vehicles in each group
	group_conf.num_veh.resize(1); group_conf.num_veh << 15;
	// Initial positions/velocities for the vehicle groups
	group_conf.veh_init.x.resize(1); group_conf.veh_init.x << 0;
	group_conf.veh_init.y.resize(1); group_conf.veh_init.y << 0;
	group_conf.veh_init.z.resize(1); group_conf.veh_init.z << 0;
	group_conf.veh_init.u.resize(1); group_conf.veh_init.u << 0;
	group_conf.veh_init.v.resize(1); group_conf.veh_init.v << 0;
	group_conf.veh_init.w.resize(1); group_conf.veh_init.w << 0;
	// Smoothing width for each group
	group_conf.veh_h.resize(1);	group_conf.veh_h << 2;
	// Limits for speed and acceleration
	group_conf.veh_limits.vmin.resize(1); group_conf.veh_limits.vmin << 3;
	group_conf.veh_limits.vmax.resize(1); group_conf.veh_limits.vmax << 6;
	group_conf.veh_limits.turning_radius.resize(1); group_conf.veh_limits.turning_radius << 0.1;
	// Total number of obstacle particles
	group_conf.num_obs = 5;
	// Size of obstacles particles * 0.5
	group_conf.obs_h.resize(6); group_conf.obs_h << 2,2,2,2,2,2;
	// Positions for the obstacles
	group_conf.obs_init.x.resize(5); group_conf.obs_init.x << 7, 12, 16,  9, 22;
	group_conf.obs_init.y.resize(5); group_conf.obs_init.y << 0,  4,  2, -2, -6;
	group_conf.obs_init.z.resize(5); group_conf.obs_init.z << 0,  0,  0,  0,  0;
	// Total number of reduced density particles
	group_conf.num_rd = 0;
	// The group that each reduced density particle belongs to.
	// Group number corresponds to array index for num_veh. -1 means not active. // NOTE: octave uses 0 for non-active!
	group_conf.rd_group.resize(1); group_conf.rd_group << 0;
	// Initial positions/velocities for reduced density particles
	group_conf.rd_init.x.resize(1); group_conf.rd_init.x << 0;
	group_conf.rd_init.y.resize(1); group_conf.rd_init.y << 0;
	group_conf.rd_init.z.resize(1); group_conf.rd_init.z << 0;
	group_conf.rd_init.u.resize(1); group_conf.rd_init.u << 0;
	group_conf.rd_init.v.resize(1); group_conf.rd_init.v << 0;
	group_conf.rd_init.w.resize(1); group_conf.rd_init.w << 0;
	// Smoothing width for reduced density particle
	group_conf.rd_h.resize(1); group_conf.rd_h << 30;
	// Total number of loiter circles
	group_conf.num_loiter = 1;
	// The group that each loiter circle belongs to.
	// Group number corresponds to array index for num_veh. -1 means not active. // NOTE: octave uses 0 for non-active!
	group_conf.loiter_group.resize(1); group_conf.loiter_group << 0;
}

/// Resize property matrices
/**
Resize prop matrices to match number of particles (#npart).

@return n/a
*/
void sph_sim::resize_prop() {
	int N = group_conf.num_veh.sum() + group_conf.num_obs + group_conf.num_rd;
	prop.vmin.resize(N, 1);
	prop.vmax.resize(N, 1);
	prop.turning_radius.resize(N, 1);
	prop.amax.resize(N, 1);
	prop.h.resize(N, 1);
	prop.m.resize(N, 1);
	prop.mu.resize(N, 1);
	prop.K.resize(N, 1);
	prop.group.resize(N, 1);
	prop.particle_type.resize(N, 1);
}

/// Initialize SPH properties
/**
Initialize #prop given the parameters specified by #group_conf and #param.

@return n/a
*/
void sph_sim::init_prop() {
	int N = 0;
	resize_prop();

	// Initialize vehicles
	for(int i = 0; i < group_conf.num_veh.size(); ++i) {
		for(int j = 0; j < group_conf.num_veh(i); ++j) {
			// Motion constraints
			prop.vmin(N,0) = group_conf.veh_limits.vmin(i);
			prop.vmax(N,0) = group_conf.veh_limits.vmax(i);
			prop.turning_radius(N,0) = group_conf.veh_limits.turning_radius(i);
			prop.amax(N,0) = pow(prop.vmax(N),2.0) / prop.turning_radius(N);

			// Smoothing width
			prop.h(N,0) = group_conf.veh_h(i);

			// Mass
			prop.m(N,0) = rho0 / kernel(0, prop.h(N), 2);

			// Kernel values at 0 and h
			double KER0 = kernel(0, prop.h(N), 2);
			double KERh = kernel(prop.h(N), prop.h(N), 2);

			// Kernel gradient at h
			double KERG = kernel_grad(prop.h(N), prop.h(N), 2);

			// Pressure force ~ K*Fp
			double Fp = rho0 * KERh / KER0;
			// Viscous force ~ mu*Fmu
			double Fmu = 2.0 * prop.vmax(N) * KER0 * KERG / ( rho0 * prop.h(N) * pow((KER0+KERh),2.0) );

			// Force coeficients found by solving:
			//	amax = vmax^2/turning_radius = K*Fp + mu*Fmu
			//	Re = (K*Fp)/(mu*Fmu)
			//
			// This enforces the desired Reynolds number at r_ij=h
			// and limits the acceleration magnitude to be on the
			// order of amax.
			prop.mu(N,0) = prop.vmax(N) / prop.turning_radius(N) / Fmu / (1.0 + param.Re);
			prop.K(N,0) = param.accel.veh * param.Re * prop.mu(N) * Fmu / Fp;

			// Group number and particle type
			prop.group(N,0) = i;
			prop.particle_type(N,0) = particle_type_enum::veh;	// NOTE: may produce oob error?
			N++;
		}
	}

	nveh = N;

	// Obstacles
	for(int i = 0; i < group_conf.num_obs; ++i) {
		// Motion constraints
		prop.vmin(N,0) = 0;
		prop.vmax(N,0) = 0;
		prop.turning_radius(N,0) = 0;
		prop.amax(N,0) = 0;

		// Smoothing width
		prop.h(N,0) = group_conf.obs_h(i);

		// Mass
		prop.m(N,0) = 2.0 * rho0 / kernel(0, prop.h(N), 2);

		// Kernel values at 0 and h
		double KER0 = kernel(0, prop.h(N), 2);
		double KERh = kernel(prop.h(N), prop.h(N), 2);

		// Kernel gradient at h
		//double KERG = kernel_grad(prop.h(N), prop.h(N), 2); // NOTE: not used?

		// Force coeficients:
		prop.mu(N,0) = 0;
		prop.K(N,0) = param.accel.obs * prop.amax.maxCoeff() * KER0 / (rho0 * KERh);

		// Group number and particle type
		prop.group(N,0) = -1;
		prop.particle_type(N,0) = particle_type_enum::obs;
		N++;
	}

	nobs = group_conf.num_obs;

	// Reduced density particles
	for(int i = 0; i < group_conf.num_rd; ++i) {
		// Motion constraints
		prop.vmin(N,0) = 0;
		prop.vmax(N,0) = 0;
		prop.turning_radius(N,0) = 0;
		prop.amax(N,0) = 0;

		// Smoothing width
		prop.h(N,0) = group_conf.rd_h(i);

		// Mass
		prop.m(N,0) = rho0 / kernel(0, prop.h(N), 1) * 1e-8;

		// Force coeficients:
		// No viscosity for attractors
		prop.mu(N,0) = 0;
		Eigen::MatrixXi I = find( ( prop.group(Eigen::seq(0,nveh-1),0).array() == group_conf.rd_group(i) ).cast<double>() );
		double max_amax = I.unaryExpr([&](int x) {
			return prop.amax(x);
		}).maxCoeff();
		prop.K(N,0) = -1.0 * param.accel.rd * max_amax * kernel(0, prop.h(N),1) / kernel_grad(prop.h(N),prop.h(N), 1);

		// Group number and particle type
		prop.group(N,0) = group_conf.rd_group(i);
		prop.particle_type(N,0) = particle_type_enum::rd;
		N++;
	}
	nrd = group_conf.num_rd;
	// Total number of SPH particles
	npart = N;
}

/// Compute h_ij matrix
/**
Set smoothing width h for each particle.

@return n/a
*/
void sph_sim::compute_hij() {
	Eigen::MatrixXd hi = prop.h(Eigen::seq(0,nveh-1),0) * Eigen::MatrixXd::Ones(1,nveh);
	Eigen::MatrixXd hj = hi.transpose();

	// Vehicles: hij=max(hi,hj)
	prop.hij = ( hi.array() > hj.array() ).select(hi,hj);

	// Obstacles: hij=h_obs
	int I1 = nveh; // NOTE: +1 for 1 based indexing? Maybe not necessary. Consider I2 as well
	int I2 = nveh + group_conf.num_obs - 1;
	prop.hij = append_down(prop.hij, prop.h(Eigen::seq(I1,I2),0) * Eigen::MatrixXd::Ones(1,prop.hij.cols()));
	prop.hij = append_right(prop.hij, Eigen::MatrixXd::Ones(prop.hij.rows(),1) * prop.h(Eigen::seq(I1,I2),0).transpose());

	// Reduced density particles: hij=h_rd
	I1 = nveh + group_conf.num_obs;
	I2 = nveh + group_conf.num_obs + group_conf.num_rd - 1;
	prop.hij = append_down(prop.hij, prop.h(Eigen::seq(I1,I2),0) * Eigen::MatrixXd::Ones(1,prop.hij.cols()));
	prop.hij = append_right(prop.hij, Eigen::MatrixXd::Ones(prop.hij.rows(),1) * prop.h(Eigen::seq(I1,I2),0).transpose());
}


/// Create a matrix that tells which kernel to use for particle interactions
/**
1 is for vehicle-reduced density particle interactions

2 is for all others

@return n/a
*/
void sph_sim::kernel_type() {
	int N = max(prop.m.rows(),prop.m.cols());

	Eigen::MatrixXd ki = prop.particle_type * Eigen::MatrixXd::Ones(1,N);
	Eigen::MatrixXd kj = ki.transpose();

	// Set everything to type 2 initially
	prop.kernel_type = 2*Eigen::MatrixXd::Ones(N,N);

	// Set vehicle-reduced density particle interactions to type 1...
	Eigen::MatrixXd lhs = ( ki.array() == (int)particle_type_enum::veh*Eigen::MatrixXd::Ones(N,N).array() ).cast<double>();
	Eigen::MatrixXd rhs = ( kj.array() == (int)particle_type_enum::rd*Eigen::MatrixXd::Ones(N,N).array() ).cast<double>();
	Eigen::MatrixXi I = find( (lhs.array() != 0 && rhs.array() != 0).cast<double>() );
	for(int i = 0; i < I.size(); ++i) {
		prop.kernel_type(I(i)) = 1;
	}
	// ...and vice-versa
	lhs = ( kj.array() == (int)particle_type_enum::veh*Eigen::MatrixXd::Ones(N,N).array() ).cast<double>();
	rhs = ( ki.array() == (int)particle_type_enum::rd*Eigen::MatrixXd::Ones(N,N).array() ).cast<double>();
	I = find( (lhs.array() != 0 && rhs.array() != 0).cast<double>() );
	for(int i = 0; i < I.size(); ++i) {
		prop.kernel_type(I(i)) = 1;
	}
}

/// Set the initial SPH states (positions and velocities) for all particles
/**
@return n/a
*/
void sph_sim::init_states() {
	if(param.ndim == 2) {
		// 2D initialization
		init2d();
	} else {
		// 3D initialization
		init3d();
	}

	// Obstacles
	for(int i = 0; i < group_conf.num_obs; ++i) {
		Eigen::MatrixXd app(1,6);
		app << group_conf.obs_init.x(i), group_conf.obs_init.y(i), group_conf.obs_init.z(i), 0, 0, 0;	// x,y,z, u,v,w
		states = append_down(states, app);
	}

	// Reduced density particles
	states = append_down(states, Eigen::MatrixXd::Zero(group_conf.num_rd,6));

}

/// Initialize particle positions for 2D simulation
/**
Uses a 2D hexagonal (close packed) lattice

@return n/a
*/
void sph_sim::init2d() {
	// Basis vectors for a hexagonal lattice in 2D
	Eigen::Vector3d v1(1, 0, 0);
	Eigen::Vector3d v2(cos(M_PI/3.0), sin(M_PI/3.0), 0);

	int N = 90;
	// Row 1 in x-direction
	Eigen::MatrixXd r(3,N);
	r << vseq(0,N-1), Eigen::MatrixXd::Zero(2,N);
	int n = r.cols();
	for(int i = 1; i < N; ++i) {	// NOTE: Check that this is
		Eigen::MatrixXd app(r.rows(), n);
		if(i % 2 == 0) {
			app << r(0,Eigen::seq(Eigen::last-n+1,Eigen::last)).array()+v2(0),
				   r(1,Eigen::seq(Eigen::last-n+1,Eigen::last)).array()+v2(1),
				   r(2,Eigen::seq(Eigen::last-n+1,Eigen::last)).array()+v2(2);
		} else {
			app << r(0,Eigen::seq(Eigen::last-n+1,Eigen::last)).array()+v2(0)-v1(0),
				   r(1,Eigen::seq(Eigen::last-n+1,Eigen::last)).array()+v2(1)-v1(1),
				   r(2,Eigen::seq(Eigen::last-n+1,Eigen::last)).array()+v2(2)-v1(2);
		}
		r = append_right(r, app);
	}

	// Randomize slightly to avoid singularities
	r(Eigen::seq(0,1),Eigen::all) = r(Eigen::seq(0,1),Eigen::all).array() + Eigen::MatrixXd::Random(2,r.cols()).array()/2*1e-8;

	// Shift to origin
	Eigen::MatrixXd ave = r.rowwise().mean();
	r.row(0) = r.row(0).array()-ave(0);
	r.row(1) = r.row(1).array()-ave(1);
	r.row(2).setZero();

	// Sort r by distance from the origin
	// r is a 3xN matrix, with x y z for the rows
	Eigen::MatrixXd d = r.array().pow(2).colwise().sum().sqrt();
	Eigen::RowVectorXi I;
	I.setLinSpaced(d.size(), 0, d.size()-1);
	sort(I.data(), I.data()+I.size(), [&](int i, int j) { return d(i) < d(j); } );
	Eigen::MatrixXd temp = r;
	r = temp(Eigen::all, I);


	// Shift to origin
	ave = r.col(0);
	r.row(0) = r.row(0).array()-ave(0);
	r.row(1) = r.row(1).array()-ave(1);
	r.row(2).setZero();

	// Sort r by distance from the origin
	d = r.array().pow(2).colwise().sum().sqrt();
	I.setLinSpaced(d.size(), 0, d.size()-1);
	sort(I.data(), I.data()+I.size(), [&](int i, int j) { return d(i) < d(j); } );
	temp = r;
	r = temp(Eigen::all, I);

	// Shift to origin
	r.row(0) = r.row(0).array()-r(0,0);
	r.row(1) = r.row(1).array()-r(0,1);
	r.row(2).setZero();

	r.transposeInPlace();
	states.resize(0,6);

	// Set the initial positions
	for(int i = 0; i < group_conf.veh_init.x.size(); ++i) {
		int n;
		bool one_loop;
		if(group_conf.veh_init.x.size() < group_conf.num_veh.size()) {
			n = group_conf.num_veh.sum();
			one_loop = true;
		} else {
			n = group_conf.num_veh(i);
			one_loop = false;
		}

		double xmid = group_conf.veh_init.x(i);
		double ymid = group_conf.veh_init.y(i);
		double zmid = 0;

		double vx = group_conf.veh_init.u(i);
		double vy = group_conf.veh_init.v(i);
		double vz = 0;

		Eigen::MatrixXd rtmp = Eigen::MatrixXd::Zero(r.rows(), 6);
		rtmp.topLeftCorner(r.rows(),3) = r(Eigen::all,Eigen::seq(0,2)) * 2 * group_conf.veh_h(i);
		rtmp.col(3).setConstant(vx);
		rtmp.col(4).setConstant(vy);
		rtmp.col(5).setZero();
		
		Eigen::MatrixXd rn = rtmp(Eigen::seqN(0,n),Eigen::all);

		double xave = rn.col(0).mean();
		double yave = rn.col(1).mean();
		double zave = 0;

		rn.col(0) = rn.col(0).array()-xave+xmid;
		rn.col(1) = rn.col(1).array()-yave+ymid;
		rn.col(2).setZero();

		
		states = append_down(states, rn);

		if(one_loop) {
			break;
		}
	}
}

/// Initialize particle positions for 3D simulation
/**
Uses a 3D hexagonal (close packed) lattice

@return n/a
*/
void sph_sim::init3d() {
	// Basis vectors for a hexagonal lattice in 3D
	//Vector3d v1(1, 0, 0);
	//Vector3d v2(cos(M_PI/3.0), sin(M_PI/3.0), 0);
	// implement later
}

/// Compute the right hand side of the SPH momentum equation
/**

\f$
\begin{align*}
\frac{D\vec{v_i}^*}{Dt^*}=-\frac{1}{M^2}\sum_{j}{ m_j^* \left ( \frac{P_i^*}{\rho_i^{*2}} + \frac{P_j^*}{\rho_j^{*2}} \right ) \frac{\partial W_{ij}^*}{\partial \vec{x_i}^*} } - \frac{1}{Re} \sum_{j}{ m_j^* \frac{\vec{\Pi_{ij}}^*}{r_{ij}^*} \frac{\partial W_{ij}^*}{\partial r_{ij}^*} }
\end{align*}
\f$

@return sph		Matrix of forces
*/
Eigen::MatrixXd sph_sim::sph_rhs() {
	// Compute the interparticle distance and vectors
	Eigen::MatrixXd dij;
	Matrix3D rij;
	Matrix3D unit_ij;
	tie(dij, rij, unit_ij) = sph_compute_dij();

	// Compute the mask function Mask(i,j) = 1 if particles interact
	Eigen::MatrixXd Mask;
	Eigen::MatrixXi MaskI;
	tie(Mask, MaskI) = sph_compute_mask(dij);

	// Compute density
	Eigen::MatrixXd rho = sph_compute_density(dij, Mask, MaskI);	// NOTE: eq (3)

	// Remove diagonal elements from MaskI
	Eigen::MatrixXd dm = Eigen::MatrixXd::Zero(npart,npart);
	dm.diagonal() = Eigen::VectorXd::Ones(npart);
	MaskI = find( ( (Mask - dm).array() == 1 ).cast<double>() );	// NOTE: can't we just set it to 0?
	
	// Compute gradW
	Eigen::MatrixXd gradW = Eigen::MatrixXd::Zero(npart, npart);
	for(int i = 0; i < MaskI.size(); ++i) {
		gradW(MaskI(i)) = kernel_grad(dij(MaskI(i)), prop.hij(MaskI(i)), prop.kernel_type(MaskI(i)));
	}

	// Compute pressure
	Eigen::MatrixXd P = sph_compute_pressure(rho);		// NOTE: eq (5)/(8)
	Eigen::MatrixXd P_term = (P.array() / rho.array().pow(2)).matrix() * prop.m.transpose() + Eigen::MatrixXd::Ones(npart,1) * (P.array() * prop.m.array() / rho.array().pow(2)).transpose().matrix();
	// Magnitude of the pressure force
	P_term = P_term.array() * gradW.array();

	// Viscosity
	Eigen::MatrixXd Pi_term = sph_compute_pi(rho, dij, rij, unit_ij, gradW, Mask, MaskI);

	Eigen::MatrixXd DvDt(unit_ij.z0.rows(), 3);
	DvDt <<	(P_term.array() * unit_ij.z0.array() * -1).rowwise().sum(),
			(P_term.array() * unit_ij.z1.array() * -1).rowwise().sum(),
			(P_term.array() * unit_ij.z2.array() * -1).rowwise().sum();
	DvDt = DvDt.array() + Pi_term.array()/param.Re;

	// External forcing
	Eigen::MatrixXd Fx, Fy, Fz;
	tie(Fx,Fy,Fz) = external_force();

	// Sum of forces
	double max_amax = Eigen::MatrixXd::NullaryExpr(nveh, 1, [&](Eigen::Index i) {
		return prop.amax(i);
	}).maxCoeff();
	DvDt(Eigen::all,0) = param.gain.sph * DvDt(Eigen::all,0).array() + param.gain.ext * max_amax * Fx.array() - param.gain.drag * states(Eigen::all,3).array();
	DvDt(Eigen::all,1) = param.gain.sph * DvDt(Eigen::all,1).array() + param.gain.ext * max_amax * Fy.array() - param.gain.drag * states(Eigen::all,4).array();
	DvDt(Eigen::all,2) = param.gain.sph * DvDt(Eigen::all,2).array() + param.gain.ext * max_amax * Fz.array() - param.gain.drag * states(Eigen::all,5).array();

	Eigen::MatrixXd rhs = sph_compute_rates(DvDt);

	if(param.ndim == 2) {
		rhs(Eigen::all,2).array() = 0;
		rhs(Eigen::all,5).array() = 0;
	}

	return rhs;
}

/// Compute the distance, vector, and unit vector between particles i and j
/**

@return 	Returns a <Eigen::MatrixXd, Matrix3D, Matrix3D> tuple, holding the distance, vector, and unit vector matrices, respectively.
*/
tuple<Eigen::MatrixXd, Matrix3D, Matrix3D> sph_sim::sph_compute_dij() {
	// Create distance matrix for dij(i,j) = distance between particles i and j
	Eigen::MatrixXd dx = states(Eigen::all,0) * Eigen::MatrixXd::Ones(1,npart);
	Eigen::MatrixXd temp = dx.transpose();
	dx = dx-temp;
	Eigen::MatrixXd dy = states(Eigen::all,1) * Eigen::MatrixXd::Ones(1,npart);
	temp = dy.transpose();
	dy = dy-temp;
	Eigen::MatrixXd dz = states(Eigen::all,2) * Eigen::MatrixXd::Ones(1,npart);
	temp = dz.transpose();
	dz = dz-temp;

	Eigen::MatrixXd dij = (dx.array().pow(2) + dy.array().pow(2) + dz.array().pow(2)).sqrt();

	Matrix3D rij(dx,dy,dz);
	Matrix3D unit_ij(dx.array()/dij.array(), dy.array()/dij.array(), dz.array()/dij.array());
	// Get rid of Infs generated by divisions by 0
	// NOTE: Check this
	unit_ij.z0 = (unit_ij.z0.array().isNaN()).select(0, unit_ij.z0);
	unit_ij.z1 = (unit_ij.z1.array().isNaN()).select(0, unit_ij.z1);
	unit_ij.z2 = (unit_ij.z2.array().isNaN()).select(0, unit_ij.z2);
	return make_tuple(dij, rij, unit_ij);
}

/// Compute the masking function and the indices that are non-zero
/**
			 { 0 if dij>2*hij
	M(i,j) = { 0 if particle j is not a vehicle and group(i)~=group(j)
			 { 0 if particle i is not a vehicle (except M(i,i)=1)
			 { 1 else
@return		Returns a <Eigen::MatrixXd, Eigen::MatrixXi> tuple, holding the mask and a column vector of indices of nonzero elements of the mask, respectively.
*/
tuple<Eigen::MatrixXd, Eigen::MatrixXi> sph_sim::sph_compute_mask(const Eigen::MatrixXd& dij) {
	// NOTE: This function uses sparse matrices, but let's ignore that for now
	// Kernel is non-zero (i.e., dij < 2*hij)
	Eigen::MatrixXd M = (dij(Eigen::seqN(0,nveh),Eigen::all).array() < 2*prop.hij(Eigen::seqN(0,nveh),Eigen::all).array()).cast<double>();
	
	// Obstacle or reduced density particle
	// NOTE: check that I'm understanding this correctly as an append down
	M = append_down(M, Eigen::MatrixXd::Zero(npart-nveh,dij.cols()));			// M(i,:) = 0
	int n = nobs + nrd;
	M(Eigen::seq(nveh,Eigen::last),Eigen::seq(nveh,Eigen::last)).diagonal() = Eigen::VectorXd::Ones(n);	// M(i,i) = 1

	// Reduced density particles
	for(int i = 0; i < nrd; ++i) {
		int I1 = nveh + nobs + i;
		Eigen::MatrixXi I2 = find( ( prop.group.array() != prop.group(I1) ).cast<double>() );
		M(I2.reshaped(),I1).array() = 0; // NOTE: check this
	}

	// Indices of nonzeros
	Eigen::MatrixXi I = find( ( M.array() != 0 ).cast<double>() );
	
	return make_tuple(M,I);
}


/// Compute the particle density equation
/**

\f$
\begin{align*}
\rho_i=\sum_j{ W ( \vec{r_{ij}} , h ) m_j }
\end{align*}
\f$

@return Matrix of rho values 
*/
Eigen::MatrixXd sph_sim::sph_compute_density(const Eigen::MatrixXd& dij, const Eigen::MatrixXd& Mask, const Eigen::MatrixXi& MaskI) {
	// Reshape mask vector into matrix
	Eigen::MatrixXd mj = Mask.array() * (Eigen::MatrixXd::Ones(npart,1) * prop.m.transpose()).array();
	
	// NOTE: Another sparse matrix
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(npart,npart);
	for(int i = 0; i < MaskI.size(); ++i) {
		K(MaskI(i)) = kernel(dij(MaskI(i)), prop.hij(MaskI(i)), prop.kernel_type(MaskI(i)));
	}

	Eigen::MatrixXd rho = ( mj.array()*K.array() ).rowwise().sum();

	// Reduced density particles have fixed density that does not consider the proximity of other particles
	Eigen::MatrixXi I = vseq(nveh+nobs, npart-1).cast<int>();
	for(int i = 0; i < I.size(); ++i) {
		rho(I(i)) = prop.m(I(i)) * kernel(0, prop.h(I(i)), 2);
	}
	
	return rho;
}

/// Equation of state to compute the pressure
/**

\f$
\begin{align*}
P=B\left(\frac{\rho}{\rho_0}-1\right)
\end{align*}
\f$

Where B is the bulk modulus.

This is the first term of the DvDt equation.

@param	rho		Particle density

@return P		Pressure
*/
Eigen::MatrixXd sph_sim::sph_compute_pressure(const Eigen::MatrixXd& rho) {
	Eigen::MatrixXd P = prop.K.array() * rho.array() * (rho.array() / rho0 - 1);

	return P;
}

///
/**
Equation to compute the viscous force.

\f$
\begin{align*}
\left\|F_\mu\right\|=\mu \left| \frac{2mv_{max}}{\rho^2 \vec{r_{i j}}} \frac{dW_1(\vec{r_{i j}},h_{ij})}{d\left\|\vec{r_{i j}}\right|} \right|
\end{align*}
\f$

This is the second term of the DvDt equation.

@param	rho		Particle density
@param	dij		Distance to particles
@param	rij		Distance to particles, dx,dy,dz
@param	unit_ij	Unit vector version of rij
@param	gradW	Kernel gradient
@param	Mask	See sph_compute_mask()
@param	MaskI	See sph_compute_mask(). Diagonal has been set to 0.

@return Pi		Viscous force
*/
Eigen::MatrixXd sph_sim::sph_compute_pi(const Eigen::MatrixXd& rho, const Eigen::MatrixXd& dij, const Matrix3D& rij, const Matrix3D& unit_ij,
								const Eigen::MatrixXd& gradW, const Eigen::MatrixXd& Mask, const Eigen::MatrixXi& MaskI) {
	Eigen::MatrixXd tmp = ( rho.array().pow(-1).matrix() * (2 * prop.m.array() / rho.array()).transpose().matrix() ).transpose().array() * gradW.array();
	for(int i = 0; i < MaskI.size(); ++i) {
		tmp(MaskI(i)) = tmp(MaskI(i))/dij(MaskI(i));
	}

	Matrix3D vji( Eigen::MatrixXd::Ones(npart,1) * states(Eigen::all,3).transpose() - states(Eigen::all,3) * Eigen::MatrixXd::Ones(1,npart),
				  Eigen::MatrixXd::Ones(npart,1) * states(Eigen::all,4).transpose() - states(Eigen::all,4) * Eigen::MatrixXd::Ones(1,npart),
				  Eigen::MatrixXd::Ones(npart,1) * states(Eigen::all,5).transpose() - states(Eigen::all,5) * Eigen::MatrixXd::Ones(1,npart) );
	
	// No viscosity for reduced density particles or obstacles
	vji.z0(Eigen::all,Eigen::seq(nveh,Eigen::last)) = Eigen::MatrixXd::Zero(vji.z0.rows(), vji.z0.cols()-nveh);	// NOTE: Check cols is correct on Zero matrix
	vji.z1(Eigen::all,Eigen::seq(nveh,Eigen::last)) = Eigen::MatrixXd::Zero(vji.z1.rows(), vji.z1.cols()-nveh);
	vji.z2(Eigen::all,Eigen::seq(nveh,Eigen::last)) = Eigen::MatrixXd::Zero(vji.z2.rows(), vji.z2.cols()-nveh);

	Eigen::MatrixXd Pi(vji.z0.rows(), 3);
	Pi <<	(tmp.array() * vji.z0.array() * -1).rowwise().sum(),
			(tmp.array() * vji.z1.array() * -1).rowwise().sum(),
			(tmp.array() * vji.z2.array() * -1).rowwise().sum();
	return Pi;
}

/// Compute the external force on vehicles to drive them toward a loiter circle
/**
@return		Returns a <Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd> tuple containing the Fx,Fy,Fz forces, respectively.
*/
tuple<Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd> sph_sim::external_force() {
	Eigen::MatrixXd Fx = Eigen::MatrixXd::Zero(states.rows(),1);
	Eigen::MatrixXd Fy = Fx;
	Eigen::MatrixXd Fz = Fx;
	for(int i = 0; i < group_conf.num_loiter; ++i) {
		int group_num = group_conf.loiter_group(i);
		Eigen::MatrixXi II = find( ( prop.group.array() == group_num ).cast<double>() );
		
		// Loiter circle
		if(lR(i) > 0) {
			// Width of the "flat spot" in the potential field, controls the width of the loiter circle track
			double width = lR(i)/4;

			// Shift the center of the loiter circle
			Eigen::MatrixXd x = states(II.reshaped(),0).array() - lx(i,0);
			Eigen::MatrixXd y = states(II.reshaped(),1).array() - lx(i,1);

			// Attraction component
			Eigen::MatrixXd d = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			d = (d.array() - lR(i)) / width;
			Eigen::MatrixXd mag = ( d.array().tanh() + d.array() / d.array().cosh().pow(2) ) * -1;

			Eigen::MatrixXd rr = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			Eigen::MatrixXd F1x = mag.array() * x.array() / rr.array();
			Eigen::MatrixXd F1y = mag.array() * y.array() / rr.array();

			// Circulation component
			Eigen::MatrixXd theta = Eigen::MatrixXd::NullaryExpr(y.rows(), y.cols(), [&](Eigen::Index i) { 
				return atan2(y(i),x(i));
			});
			Eigen::MatrixXd F2x = ( exp(2) * ( rr.array() / lR(i) ).pow(2) * exp(-2 * rr.array() / lR(i)) ) * sin(theta.array()) * -1;
			Eigen::MatrixXd F2y = ( exp(2) * ( rr.array() / lR(i) ).pow(2) * exp(-2 * rr.array() / lR(i)) ) * cos(theta.array());

			// Total force
			double w = 1.0;
			for(int i = 0; i < II.size(); ++i) {
				Fx(II(i)) = w*F1x(i) + (2-w)*F2x(i);
				Fy(II(i)) = w*F1y(i) + (2-w)*F2y(i);
			}
		} else {
			// Simple attractor (no circulation force)

			// Shift the center of the loiter circle
			double width = sqrt(
								pow(
									II.unaryExpr( [&](int x) { return prop.h(x); } ).mean(),
									2
								)
								* II.size()
							) / 2.0;
			Eigen::MatrixXd x = states(II.reshaped(),0) - lx(II.reshaped(),0);
			Eigen::MatrixXd y = states(II.reshaped(),1) - lx(II.reshaped(),1);

			// Attraction component
			Eigen::MatrixXd d = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			d = d.array() / width;
			Eigen::MatrixXd mag = ( d.array().tanh() + d.array() / d.array().cosh().pow(2) ) * -1;

			Eigen::MatrixXd rr = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			for(int i = 0; i < II.size(); ++i) {
				Fx(II(i)) = mag(i) * x(i) / rr(i);
				Fy(II(i)) = mag(i) * y(i) / rr(i);
			}
		}
	}
	Fx = Fx.array().isNaN().select(0, Fx);
	Fy = Fy.array().isNaN().select(0, Fy);

	return make_tuple(Fx,Fy,Fz);
}

/// Compute the rate of change of the velocity and accelerations, while applying vehicle constraints
/**
@param	DvDt	Matrix of the SPH momentum equation
@return	rates	Matrix of rates of change in velocity and accelerations.
*/
Eigen::MatrixXd sph_sim::sph_compute_rates(const Eigen::MatrixXd& DvDt) {
	// Break acceleration into 2 components, normal and tangential:
	Eigen::MatrixXd v = states(Eigen::all,Eigen::seq(3,5));
	Eigen::MatrixXd vmag = v.array().pow(2).rowwise().sum().sqrt();
	// Unit vector in the v direction
	Eigen::MatrixXd vhat = v.array() / (vmag * Eigen::MatrixXd::Ones(1,3) ).array();
	Eigen::MatrixXi I = find( ( vmag.array() == 0 ).cast<double>() );
	vhat(I.reshaped(),0).array() = 1;
	vhat(I.reshaped(),1).array() = 0;
	vhat(I.reshaped(),2).array() = 0;

	// Acceleration in the normal and tangent direction
	Eigen::MatrixXd a_tan = ( ( DvDt.array() * vhat.array() ).rowwise().sum().matrix() * Eigen::MatrixXd::Ones(1,3) ).array() * vhat.array();
	Eigen::MatrixXd a_norm = DvDt - a_tan;
	Eigen::MatrixXd a_tan_mag = a_tan.array().pow(2).rowwise().sum().sqrt();

	// Limit acceleration
	I = find( ( a_tan_mag.array() > prop.amax.array() ).cast<double>() );
	// NOTE: Check all these. Also check order of operations. ./ then *, or * then ./ ?
	if(I.size() != 0 ) {
		a_tan(I.reshaped(),Eigen::all) = a_tan(I.reshaped(),Eigen::all).array() / ( a_tan_mag(I.reshaped(),0) * Eigen::MatrixXd::Ones(1,3) ).array() * ( prop.amax(I.reshaped(), 0) * Eigen::MatrixXd::Ones(1,3) ).array();
	}

	// Limit speed
	I = find( ( vmag.array() > prop.vmax.array() ).cast<double>() );
	if(I.size() != 0 ) {
		a_tan(I.reshaped(), Eigen::all) = (prop.amax(I.reshaped(),0) * Eigen::MatrixXd::Ones(1,3)).array() * vhat(I.reshaped(), Eigen::all).array();
	}
	I = find( ( vmag.array() < prop.vmin.array() ).cast<double>() );
	if(I.size() != 0 ) {
		a_tan(I.reshaped(), Eigen::all) = (prop.amax(I.reshaped(),0) * Eigen::MatrixXd::Ones(1,3)).array() * vhat(I.reshaped(), Eigen::all).array();
	}

	// Limit turning radius
	Eigen::MatrixXd a_norm_mag = a_norm.array().pow(2).rowwise().sum().sqrt();
	I = find( ( a_norm_mag.array() > vmag.array().pow(2) / prop.turning_radius.array() ).cast<double>() );
	
	if(I.size() != 0 ) {
		Eigen::MatrixXd temp = vmag(I.reshaped(),0).array().pow(2).array() / prop.turning_radius(I.reshaped(),0).array();
		a_norm(I.reshaped(), Eigen::all) = a_norm(I.reshaped(),Eigen::all).array() / (a_norm_mag(I.reshaped(),0) * Eigen::MatrixXd::Ones(1,3)).array() * (temp*Eigen::MatrixXd::Ones(1,3)).array();
	}

	Eigen::MatrixXd rates = append_right( states(Eigen::all,Eigen::seq(3,5)) , ( a_tan(Eigen::all,Eigen::seq(0,2)) + a_norm(Eigen::all,Eigen::seq(0,2)) ) );

	return rates;
}

/// Apply velocity constraints
/**
@return	n/a
*/
void sph_sim::constrain_vel() {
	// Max velocity constrain
	Eigen::MatrixXd V = states(Eigen::all,Eigen::seq(3,5)).array().pow(2).rowwise().sum().sqrt();
	Eigen::MatrixXi I = find( (V.array() > prop.vmax.array() ).cast<double>() );
	if(I.size() != 0 ) {
		states(I.reshaped(), 3) = states(I.reshaped(), 3).array() / V(I.reshaped(),0).array() * prop.vmax(I.reshaped(),0).array();
		states(I.reshaped(), 4) = states(I.reshaped(), 4).array() / V(I.reshaped(),0).array() * prop.vmax(I.reshaped(),0).array();
		states(I.reshaped(), 5) = states(I.reshaped(), 5).array() / V(I.reshaped(),0).array() * prop.vmax(I.reshaped(),0).array();
	}

	// Min velocity constrain
	V = states(Eigen::all,Eigen::seq(3,5)).array().pow(2).rowwise().sum().sqrt();
	I = find( (V.array() < prop.vmin.array() ).cast<double>() );
	if(I.size() != 0 ) {
		states(I.reshaped(), 3) = states(I.reshaped(), 3).array() / V(I.reshaped(),0).array() * prop.vmin(I.reshaped(),0).array();
		states(I.reshaped(), 4) = states(I.reshaped(), 4).array() / V(I.reshaped(),0).array() * prop.vmin(I.reshaped(),0).array();
		states(I.reshaped(), 5) = states(I.reshaped(), 5).array() / V(I.reshaped(),0).array() * prop.vmin(I.reshaped(),0).array();
	}
}

// =============================================================
// ========================== GETTERS ==========================
// =============================================================
/// Return the current time in the SPH simulation
/**
@return	double
*/
double sph_sim::get_time() {
	return t;
}
/// Return the initial time for the SPH simulation
/**
@return	double
*/
double sph_sim::get_initial_time() {
	return t0;
}
/// Return the time step to be used for the SPH simulation
/**
@return	double
*/
double sph_sim::get_dt() {
	return param.dt;
}

/// Return a matrix containing the [x y z] positions and [u v w] velocities of all SPH particles.
/**
Each particle is stored in one row:

			 [ x0 y0 z0 u0 v0 w0 ]
	states = [ x1 y1 z1 u1 v1 w1 ]
			 [		 ...		]

@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_states() {
	return states;
}

/// Return the total number of particles in the simulation
/**
@return	int
*/
int sph_sim::get_npart() {
	return npart;
}
/// Return the number of vehicles in the simulation
/**
@return	int
*/
int sph_sim::get_nveh() {
	return nveh;
}
/// Return the number of obstacles in the simulation
/**
@return	int
*/
int sph_sim::get_nobs() {
	return nobs;
}
/// Return the number of reduced density (attractor) particles in the simulation
/**
@return	int
*/
int sph_sim::get_nrd() {
	return nrd;
}

/// Return a column vector containing the x positions of all the SPH particles
/**
@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_x() {
	return states.col(0);
}
/// Return a column vector containing the y positions of all the SPH particles
/**
@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_y() {
	return states.col(1);
}
/// Return a column vector containing the z positions of all the SPH particles
/**
@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_z() {
	return states.col(2);
}
/// Return a column vector containing the u velocities of all the SPH particles
/**
@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_u() {
	return states.col(3);
}
/// Return a column vector containing the v velocities of all the SPH particles
/**
@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_v() {
	return states.col(4);
}
/// Return a column vector containing the w velocities of all the SPH particles
/**
@return	Eigen::MatrixXd
*/
Eigen::MatrixXd sph_sim::get_w() {
	return states.col(5);
}

/// Returns a prop_struct containing all the properties of each SPH particle in the simulation
/**
@return prop_struct
*/
prop_struct sph_sim::get_prop() {
	return prop;
}

/// Returns a param_struct containing all the parameters used in the SPH simulation
/**
@return param_struct
*/
param_struct sph_sim::get_param() {
	return param;
}

/// Return a group_conf_struct containing the group configuration
/**
@return group_conf_struct
*/
group_conf_struct sph_sim::get_group_conf() {
	return group_conf;
}

// TODO
// tags: NOTE, INCOMPLETE, MODIFIED
// Contemplate the type of vseq. int or double? eh...
// Check that arithmetic stuff uses doubles and not ints (e.g., 2.0 not 2)
// Consider that I'm comparing using doubles (e.g., x == 1.0), which may cause issues?
// Can maybe add an overload for the find function so I don't have to cast<double>()