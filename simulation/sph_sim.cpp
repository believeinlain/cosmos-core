#include "sph_sim.h"

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

// Return vector of indices of nonzero elements. As row if A is a row vector, as a column vector otherwise
MatrixXi find(const MatrixXd& A) {
	MatrixXi idxs(1,A.size());
	int size = 0;
	for(Index i=0; i<A.size(); ++i) {
		if(A(i)) {
			idxs(0,size) = i;
			size++;
		}
	}
	if(size > 0) {
		idxs.conservativeResize(1,size);
		if(A.rows() > 1) {
			idxs.transposeInPlace();
		}
	} else {
		idxs.resize(0,0);
	}
	
	return idxs;
}

// Return a matrix whose elements are those indices of m specified by I. Necessarily, returned matrix will be the same shape as I
// Iow: indexing the elements of m with another matrix I
MatrixXd index(const MatrixXd& m, const MatrixXi& I) {
	MatrixXd out(I.rows(),I.cols());

	for(int i = 0; i < I.size(); ++i) {
		out(i) = m(I(i));
	}

	return out;
}

// Returns a row vector of doubles in rising sequence by 1 (eg: {0, 1, 2, 3, ... }), inclusive
MatrixXd vseq(int val0, int valn) {
	RowVectorXd rvxd;
	rvxd.setLinSpaced(valn-val0+1, val0, valn);
	return rvxd.matrix();
}

// Return sorted indices of c, in ascending order
MatrixXi sort(const MatrixXd& c) {
	MatrixXi idxs(1,c.size());
	for(int i = 0; i < c.size(); ++i) {
		idxs(i) = i;
	}
	sort(idxs.data(), idxs.data()+idxs.size(),[&](int i, int j) {
		return c(i) < c(j);
	});
	return idxs;
}

// Append another matrix app to matrix m, right or down
MatrixXd append_right(const MatrixXd& m, const MatrixXd& app) {
	int off = m.cols();
	MatrixXd out(m.rows(), off + app.cols());
	out.topLeftCorner(m.rows(), m.cols()) = m;
	out(seq(0,last),seq(off,last)) = app;
	return out;
}
MatrixXd append_down(const MatrixXd& m, const MatrixXd& app) {
	int off = m.rows();
	MatrixXd out(off + app.rows(), m.cols());
	out.topLeftCorner(m.rows(), m.cols()) = m;
	out(seq(off,last),seq(0,last)) = app;
	return out;
}



// sph_sim constructor
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
sph_sim::sph_sim(param_struct param, group_conf_struct group_conf, double t0 /*= 0*/) : param(param), group_conf(group_conf), t0(t0) {
	// Setup SPH properties
	init_prop();
	compute_hij();
	kernel_type();

	// Initialize positions and velocities
	init_states();
}

// Update or change the SPH properties to match the properties in the arguments param and group_conf
void sph_sim::sph_update_properties(const param_struct& param, const group_conf_struct& group) {
	// Define parameters
	this->param = param;

	// Define group configuration
	this->group_conf = group_conf;

	// Setup SPH properties
	init_prop();
	compute_hij();
}

// Take a single time-step forward in the simulation
void sph_sim::sph_sim_step(MatrixXd rdx, MatrixXd lx, MatrixXd lR) {
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
		states.col(0)(seq(I1,last)) = rdx.col(0);
		states.col(1)(seq(I1,last)) = rdx.col(1);
		if(param.ndim == 3) {
			states.col(2)(seq(I1,last)) = rdx.col(2);
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
	MatrixXd uvw = MatrixXd::Zero(states.rows(),3);

	// Fourth order Runge-Kutta timestep
	MatrixXd k1 = sph_rhs();
	k1(all,seq(0,2)) = k1(all,seq(0,2)) + uvw;
	
	sph_sim tmp = *this;
	tmp.states = tmp.states + tmp.param.dt/2.0*k1;
	MatrixXd k2 = tmp.sph_rhs();
	k2(all,seq(0,2)) = k2(all,seq(0,2)) + uvw;

	tmp = *this;
	tmp.states = tmp.states + tmp.param.dt/2.0*k2;
	MatrixXd k3 = tmp.sph_rhs();
	k3(all,seq(0,2)) = k3(all,seq(0,2)) + uvw;

	tmp = *this;
	tmp.states = tmp.states + tmp.param.dt*k3;
	MatrixXd k4 = tmp.sph_rhs();
	k4(all,seq(0,2)) = k4(all,seq(0,2)) + uvw;

	states = states + param.dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
	// Increment time
	t = t + param.dt;

	// Constrain the velocity
	constrain_vel();
}

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
		double KERG = kernel_grad(prop.h(N), prop.h(N), 2);

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
		// NOTE: check bounds on the seqs (seq(0,nveh))
		MatrixXi I = find( ( prop.group(seq(0,nveh-1),0).array() == group_conf.rd_group(i) ).cast<double>() );
		prop.K(N,0) = -1.0 * param.accel.rd * index(prop.amax,I).maxCoeff()
							* kernel(0, prop.h(N),1) / kernel_grad(prop.h(N),prop.h(N), 1);

		// Group number and particle type
		prop.group(N,0) = group_conf.rd_group(i);
		prop.particle_type(N,0) = particle_type_enum::rd;
		N++;
	}
	nrd = group_conf.num_rd;
	// Total number of SPH particles
	npart = N;
}

// Compute h_ij matrix
void sph_sim::compute_hij() { // NOTE: Check the prop.h(seq(I1,I2),0), same as above, was causing errors
	MatrixXd hi = prop.h(seq(0,nveh-1),0) * MatrixXd::Ones(1,nveh);
	MatrixXd hj = hi.transpose();

	// Vehicles: hij=max(hi,hj)
	prop.hij = ( hi.array() > hj.array() ).select(hi,hj);

	// Obstacles: hij=h_obs
	int I1 = nveh; // NOTE: +1 for 1 based indexing? Maybe not necessary. Consider I2 as well
	int I2 = nveh + group_conf.num_obs - 1;
	prop.hij = append_down(prop.hij, prop.h(seq(I1,I2),0) * MatrixXd::Ones(1,prop.hij.cols()));
	prop.hij = append_right(prop.hij, MatrixXd::Ones(prop.hij.rows(),1) * prop.h(seq(I1,I2),0).transpose());

	// Reduced density particles: hij=h_rd
	I1 = nveh + group_conf.num_obs;
	I2 = nveh + group_conf.num_obs + group_conf.num_rd - 1;
	prop.hij = append_down(prop.hij, prop.h(seq(I1,I2),0) * MatrixXd::Ones(1,prop.hij.cols()));
	prop.hij = append_right(prop.hij, MatrixXd::Ones(prop.hij.rows(),1) * prop.h(seq(I1,I2),0).transpose());
}

// Create a matrix kernel_type that tells which kernel to use.
// 1 is for vehicle-reduced density particle interactions,
// 2 is for all others
void sph_sim::kernel_type() {
	int N = prop.m.rows() > prop.m.cols() ? prop.m.rows() : prop.m.cols();

	MatrixXd ki = prop.particle_type * MatrixXd::Ones(1,N);
	MatrixXd kj = ki.transpose();

	prop.kernel_type = 2*MatrixXd::Ones(N,N);

	MatrixXd lhs = ( ki.array() == (int)particle_type_enum::veh*MatrixXd::Ones(N,N).array() ).cast<double>();
	MatrixXd rhs = ( kj.array() == (int)particle_type_enum::rd*MatrixXd::Ones(N,N).array() ).cast<double>();
	MatrixXi I = find( (lhs.array() != 0 && rhs.array() != 0).cast<double>() ); // NOTE: error here
	I = I.unaryExpr([&](int x) {
		prop.kernel_type(x) = 1;
		return x;
	});
	lhs = ( kj.array() == (int)particle_type_enum::veh*MatrixXd::Ones(N,N).array() ).cast<double>();
	rhs = ( ki.array() == (int)particle_type_enum::rd*MatrixXd::Ones(N,N).array() ).cast<double>();
	I = find( (lhs.array() != 0 && rhs.array() != 0).cast<double>() );
	I = I.unaryExpr([&](int x) {
		prop.kernel_type(x) = 1;
		return x;
	});
}

// Set the initial SPH states (positions and velocities) for all particles
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
		MatrixXd app(1,6);
		app << group_conf.obs_init.x(i), group_conf.obs_init.y(i), group_conf.obs_init.z(i), 0, 0, 0;	// x,y,z, u,v,w
		states = append_down(states, app);
	}

	// Reduced density particles
	states = append_down(states, MatrixXd::Zero(group_conf.num_rd,6));

}

// 2D initialization, use a 2D hexagonal (close packed) lattice
void sph_sim::init2d() {
	// Basis vectors for a hexagonal lattice in 2D
	Vector3d v1(1, 0, 0);
	Vector3d v2(cos(M_PI/3.0), sin(M_PI/3.0), 0);

	int N = 90;
	// Row 1 in x-direction
	MatrixXd r(3,N);
	r << vseq(0,N-1), MatrixXd::Zero(2,N);
	int n = r.cols();
	for(int i = 1; i < N; ++i) {	// NOTE: Check that this is
		MatrixXd app(r.rows(), n);
		if(i % 2 == 0) {
			app << r(0,seq(last-n+1,last)).array()+v2(0),
				   r(1,seq(last-n+1,last)).array()+v2(1),
				   r(2,seq(last-n+1,last)).array()+v2(2);
		} else {
			app << r(0,seq(last-n+1,last)).array()+v2(0)-v1(0),
				   r(1,seq(last-n+1,last)).array()+v2(1)-v1(1),
				   r(2,seq(last-n+1,last)).array()+v2(2)-v1(2);
		}
		r = append_right(r, app);
	}

	// Randomize slightly to avoid singularities
	r(seq(0,1),all) = r(seq(0,1),all).array() + MatrixXd::Random(2,r.cols()).array()/2*1e-8;

	// Shift to origin
	MatrixXd ave = r.rowwise().mean();
	r.row(0) = r.row(0).array()-ave(0);
	r.row(1) = r.row(1).array()-ave(1);
	r.row(2).setZero();

	// Sort r by distance from the origin
	// r is a 3xN matrix, with x y z for the rows
	MatrixXd d = r.array().pow(2).colwise().sum().sqrt();
	RowVectorXi I;
	I.setLinSpaced(d.size(), 0, d.size()-1);
	sort(I.data(), I.data()+I.size(), [&](int i, int j) { return d(i) < d(j); } );
	MatrixXd temp = r;
	r = temp(all, I);


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
	r = temp(all, I);

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

		MatrixXd rtmp = MatrixXd::Zero(r.rows(), 6);
		rtmp.topLeftCorner(r.rows(),3) = r(all,seq(0,2)) * 2 * group_conf.veh_h(i);
		rtmp.col(3).setConstant(vx);
		rtmp.col(4).setConstant(vy);
		rtmp.col(5).setZero();
		
		MatrixXd rn = rtmp(seqN(0,n),all);

		double xave = rn.col(0).mean();
		double yave = rn.col(1).mean();
		double zave = 0;

		rn.col(0) = rn.col(0).array()-xave+xmid;
		rn.col(1) = rn.col(1).array()-yave+xmid;
		rn.col(2).setZero();

		
		states = append_down(states, rn);

		if(one_loop) {
			break;
		}
	}
}

// 3D initialization, use a 3D hexagonal (close packed) lattice
void sph_sim::init3d() {
	// Basis vectors for a hexagonal lattice in 3D
	Vector3d v1(1, 0, 0);
	Vector3d v2(cos(M_PI/3.0), sin(M_PI/3.0), 0);
	// implement later
}

// Return the right hand side of the SPH momentum equation
MatrixXd sph_sim::sph_rhs() {
	// Compute the interparticle distance and vectors
	MatrixXd dij;
	Matrix3D rij;
	Matrix3D unit_ij;
	tie(dij, rij, unit_ij) = sph_compute_dij();

	// Compute the mask function Mask(i,j) = 1 if particles interact
	MatrixXd Mask;
	MatrixXi MaskI;
	tie(Mask, MaskI) = sph_compute_mask(dij);

	// Compute density
	MatrixXd rho = sph_compute_density(dij, Mask, MaskI);

	// Remove diagonal elements from MaskI
	DiagonalMatrix<double, Dynamic> dm(npart);
	dm.diagonal() = VectorXd::Ones(npart);
	MatrixXd tmp = Mask - (MatrixXd)dm;	// NOTE: Check this
	MaskI = find( ( tmp.array() == 1 ).cast<double>() );
	
	// Compute gradW
	MatrixXd gradW = MatrixXd::Zero(npart, npart);
	MaskI = MaskI.unaryExpr([&](int x) {
		gradW(x) = kernel_grad(dij(x), prop.hij(x), prop.kernel_type(x));
		return x;
	});

	// Compute pressure
	MatrixXd P = sph_compute_pressure(rho);
	// NOTE: Check this
	MatrixXd P_term = (P.array() / rho.array().pow(2)).matrix() * prop.m.transpose() + MatrixXd::Ones(npart,1) * (P.array() * prop.m.array() / rho.array().pow(2)).transpose().matrix();
	// Magnitude of the pressure force
	P_term = P_term.array() * gradW.array();

	// Viscosity
	MatrixXd Pi_term = sph_compute_pi(rho, dij, rij, unit_ij, gradW, Mask, MaskI);

	MatrixXd DvDt(unit_ij.z0.rows(), 3);
	DvDt <<	(P_term.array() * unit_ij.z0.array() * -1).rowwise().sum(),
			(P_term.array() * unit_ij.z1.array() * -1).rowwise().sum(),
			(P_term.array() * unit_ij.z2.array() * -1).rowwise().sum();
	DvDt = DvDt.array() + Pi_term.array()/param.Re;

	// External forcing
	MatrixXd Fx, Fy, Fz;
	tie(Fx,Fy,Fz) = external_force();

	// Sum of forces
	double max_amax = MatrixXd::NullaryExpr(nveh, 1, [&](Index i) {
		return prop.amax(i);
	}).maxCoeff();
	DvDt(all,0) = param.gain.sph * DvDt(all,0).array() + param.gain.ext * max_amax * Fx.array() - param.gain.drag * states(all,3).array();
	DvDt(all,1) = param.gain.sph * DvDt(all,1).array() + param.gain.ext * max_amax * Fy.array() - param.gain.drag * states(all,4).array();
	DvDt(all,2) = param.gain.sph * DvDt(all,2).array() + param.gain.ext * max_amax * Fz.array() - param.gain.drag * states(all,5).array();

	MatrixXd rhs = sph_compute_rates(DvDt);

	if(param.ndim == 2) {
		rhs(all,2).array() = 0;
		rhs(all,5).array() = 0;
	}

	return rhs;
}

// Compute the distance, vector, and unit vector between particles i and j
tuple<MatrixXd, Matrix3D, Matrix3D> sph_sim::sph_compute_dij() {
	// Create distance matrix for dij(i,j) = distance between particles i and j
	MatrixXd dx = states(all,0) * MatrixXd::Ones(1,npart);
	MatrixXd temp = dx.transpose();
	dx = dx-temp;
	MatrixXd dy = states(all,1) * MatrixXd::Ones(1,npart);
	temp = dy.transpose();
	dy = dy-temp;
	MatrixXd dz = states(all,2) * MatrixXd::Ones(1,npart);
	temp = dz.transpose();
	dz = dz-temp;

	MatrixXd dij = (dx.array().pow(2) + dy.array().pow(2) + dz.array().pow(2)).sqrt();

	Matrix3D rij(dx,dy,dz);
	Matrix3D unit_ij(dx.array()/dij.array(), dy.array()/dij.array(), dz.array()/dij.array());
	// Get rid of Infs generated by divisions by 0
	// NOTE: Check this
	unit_ij.z0 = (unit_ij.z0.array().isNaN()).select(0, unit_ij.z0);
	unit_ij.z1 = (unit_ij.z1.array().isNaN()).select(0, unit_ij.z1);
	unit_ij.z2 = (unit_ij.z2.array().isNaN()).select(0, unit_ij.z2);
	return make_tuple(dij, rij, unit_ij);
}

// Compute the masking function and the indices that are non-zero
//			{ 0 if dij>2*hij
// M(i,j) = { 0 if particle j is not a vehicle and group(i)~=group(j)
//			{ 0 if particle i is not a vehicle (except M(i,i)=1)
//			{ 1 else
tuple<MatrixXd, MatrixXi> sph_sim::sph_compute_mask(const MatrixXd& dij) {
	// NOTE: This function uses sparse matrices, but let's ignore that for now
	// Kernel is non-zero (i.e., dij < 2*hij)
	MatrixXd M = (dij(seqN(0,nveh),all).array() < 2*prop.hij(seqN(0,nveh),all).array()).cast<double>();
	
	// Obstacle or reduced density particle
	// NOTE: check that I'm understanding this correctly as an append down
	M = append_down(M, MatrixXd::Zero(npart-nveh,dij.cols()));			// M(i,:) = 0
	int n = nobs + nrd;
	M(seq(nveh,last),seq(nveh,last)).diagonal() = VectorXd::Ones(n);	// M(i,i) = 1

	// Reduced density particles
	for(int i = 0; i < nrd; ++i) {
		int I1 = nveh + nobs + i;
		MatrixXi I2 = find( ( prop.group.array() != prop.group(I1) ).cast<double>() );
		M(I2.reshaped(),I1).array() = 0; // NOTE: check this
	}

	// Indicies of nonzeros
	MatrixXi I = find( ( M.array() != 0 ).cast<double>() );
	
	return make_tuple(M,I);
}

// Mask I is a column vector of indices
MatrixXd sph_sim::sph_compute_density(const MatrixXd& dij, const MatrixXd& Mask, const MatrixXi& MaskI) {
	// Reshape mask vector into matrix
	// NOTE: crashes here after a while because prop.m reaches 40x1. Mask is 20x20, npart is 20
	MatrixXd mj = Mask.array() * (MatrixXd::Ones(npart,1) * prop.m.transpose()).array();
	
	// NOTE: Another sparse matrix
	MatrixXd K = MatrixXd::Zero(npart,npart);
	MatrixXi temp = MaskI.unaryExpr([&](int x) {
		K(x) = kernel(dij(x), prop.hij(x), prop.kernel_type(x));
		return x;
	});

	MatrixXd rho = ( mj.array()*K.array() ).rowwise().sum();

	// Reduced density particles have fixed density that does not consider the proximity of other particles
	MatrixXi I = vseq(nveh+nobs, npart-1).cast<int>();
	I = I.unaryExpr([&](int x) {
		rho(x) = prop.m(x) * kernel(0, prop.h(x), 2);
		return x;
	});
	
	return rho;
}

// Equation of state to compute the pressure
MatrixXd sph_sim::sph_compute_pressure(const MatrixXd& rho) {
	MatrixXd P = prop.K.array() * rho.array() * (rho.array() / rho0 - 1);

	return P;
}

// Compute the viscous forces
MatrixXd sph_sim::sph_compute_pi(const MatrixXd& rho, const MatrixXd& dij, const Matrix3D& rij, const Matrix3D& unit_ij,
								const MatrixXd& gradW, const MatrixXd& Mask, const MatrixXi& MaskI) {
	MatrixXd tmp = ( rho.array().pow(-1).matrix() * (2 * prop.m.array() / rho.array()).transpose().matrix() ).transpose().array() * gradW.array();
	auto temp = MaskI.unaryExpr([&](int x) {
		tmp(x) = tmp(x)/dij(x);
		return x;
	});

	Matrix3D vji( MatrixXd::Ones(npart,1) * states(all,3).transpose() - states(all,3) * MatrixXd::Ones(1,npart),
				  MatrixXd::Ones(npart,1) * states(all,4).transpose() - states(all,4) * MatrixXd::Ones(1,npart),
				  MatrixXd::Ones(npart,1) * states(all,5).transpose() - states(all,5) * MatrixXd::Ones(1,npart) );
	
	// No viscosity for reduced density particles or obstacles
	vji.z0(all,seq(nveh,last)) = MatrixXd::Zero(vji.z0.rows(), vji.z0.cols()-nveh);	// NOTE: Check cols is correct on Zero matrix
	vji.z1(all,seq(nveh,last)) = MatrixXd::Zero(vji.z1.rows(), vji.z1.cols()-nveh);
	vji.z2(all,seq(nveh,last)) = MatrixXd::Zero(vji.z2.rows(), vji.z2.cols()-nveh);

	MatrixXd Pi(vji.z0.rows(), 3);
	Pi <<	(tmp.array() * vji.z0.array() * -1).rowwise().sum(),
			(tmp.array() * vji.z1.array() * -1).rowwise().sum(),
			(tmp.array() * vji.z2.array() * -1).rowwise().sum();
	return Pi;
}

// Compute the external force on vehicles to drive them toward a loiter circle
tuple<MatrixXd,MatrixXd,MatrixXd> sph_sim::external_force() {
	MatrixXd Fx = MatrixXd::Zero(states.rows(),1);
	MatrixXd Fy = Fx;
	MatrixXd Fz = Fx;
	for(int i = 0; i < group_conf.num_loiter; ++i) {
		int group_num = group_conf.loiter_group(i);
		MatrixXi II = find( ( prop.group.array() == group_num ).cast<double>() );
		
		// Loiter circle
		if(lR(i) > 0) {
			// Width of the "flat spot" in the potential field, controls the width of the loiter circle track
			double width = lR(i)/4;

			// Shift the center of the loiter circle
			MatrixXd x = states(II.reshaped(),0).array() - lx(i,0);
			MatrixXd y = states(II.reshaped(),1).array() - lx(i,1);

			// Attraction component
			MatrixXd d = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			d = (d.array() - lR(i)) / width;
			MatrixXd mag = ( d.array().tanh() + d.array() / d.array().cosh().pow(2) ) * -1;

			MatrixXd rr = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			MatrixXd F1x = mag.array() * x.array() / rr.array();
			MatrixXd F1y = mag.array() * y.array() / rr.array();

			// Circulation component
			MatrixXd theta = MatrixXd::NullaryExpr(y.rows(), y.cols(), [&](Index i) { 
				return atan2(y(i),x(i));
			});
			MatrixXd F2x = ( exp(2) * ( rr.array() / lR(i) ).pow(2) * exp(-2 * rr.array() / lR(i)) ) * sin(theta.array()) * -1;
			MatrixXd F2y = ( exp(2) * ( rr.array() / lR(i) ).pow(2) * exp(-2 * rr.array() / lR(i)) ) * cos(theta.array());

			// Total force
			double w = 1.0;
			II = MatrixXi::NullaryExpr(II.rows(), II.cols(), [&](Index i) {
				Fx(II(i)) = w*F1x(i) + (2-w)*F2x(i);
				Fy(II(i)) = w*F1y(i) + (2-w)*F2y(i);
				return II(i);
			});
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
			MatrixXd x = states(II.reshaped(),0) - lx(II.reshaped(),0);
			MatrixXd y = states(II.reshaped(),1) - lx(II.reshaped(),1);

			// Attraction component
			MatrixXd d = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			d = d.array() / width;
			MatrixXd mag = ( d.array().tanh() + d.array() / d.array().cosh().pow(2) ) * -1;

			MatrixXd rr = ( x.array().pow(2) + y.array().pow(2) ).sqrt();
			II = MatrixXi::NullaryExpr(II.rows(), II.cols(), [&](Index i) {
				Fx(II(i)) = mag(i) * x(i) / rr(i);
				Fy(II(i)) = mag(i) * y(i) / rr(i);
				return II(i);
			});
		}
	}
	Fx = (Fx.array().isNaN()).select(0, Fx);
	Fy = (Fy.array().isNaN()).select(0, Fy);

	return make_tuple(Fx,Fy,Fz);
}

// Compute the rate of change of SPH.states, i.e., the velocity
// and accelerations, while applying vehicle constraints
MatrixXd sph_sim::sph_compute_rates(const MatrixXd& DvDt) {
	// Break acceleration into 2 components, normal and tangential:
	MatrixXd v = states(all,seq(3,5));
	MatrixXd vmag = v.array().pow(2).rowwise().sum().sqrt();
	// Unit vector in the v direction
	MatrixXd vhat = v.array() / (vmag * MatrixXd::Ones(1,3) ).array();
	MatrixXi I = find( ( vmag.array() == 0 ).cast<double>() );
	vhat(I.reshaped(),0).array() = 1;
	vhat(I.reshaped(),1).array() = 0;
	vhat(I.reshaped(),2).array() = 0;

	// Acceleration in the normal and tangent direction
	MatrixXd a_tan = ( ( DvDt.array() * vhat.array() ).rowwise().sum().matrix() * MatrixXd::Ones(1,3) ).array() * vhat.array();
	MatrixXd a_norm = DvDt - a_tan;
	MatrixXd a_tan_mag = a_tan.array().pow(2).rowwise().sum().sqrt();

	// Limit acceleration
	I = find( ( a_tan_mag.array() > prop.amax.array() ).cast<double>() );
	// NOTE: Check all these. Also check order of operations. ./ then *, or * then ./ ?
	if(I.size() != 0 ) {
		a_tan(I.reshaped(),all) = a_tan(I.reshaped(),all).array() / ( a_tan_mag(I.reshaped(),0) * MatrixXd::Ones(1,3) ).array() * ( prop.amax(I.reshaped(), 0) * MatrixXd::Ones(1,3) ).array();
	}

	// Limit speed
	I = find( ( vmag.array() > prop.vmax.array() ).cast<double>() );
	if(I.size() != 0 ) {
		a_tan(I.reshaped(), all) = (prop.amax(I.reshaped(),0) * MatrixXd::Ones(1,3)).array() * vhat(I.reshaped(), all).array();
	}
	I = find( ( vmag.array() < prop.vmin.array() ).cast<double>() );
	if(I.size() != 0 ) {
		a_tan(I.reshaped(), all) = (prop.amax(I.reshaped(),0) * MatrixXd::Ones(1,3)).array() * vhat(I.reshaped(), all).array();
	}

	// Limit turning radius
	MatrixXd a_norm_mag = a_norm.array().pow(2).rowwise().sum().sqrt();
	I = find( ( a_norm_mag.array() > vmag.array().pow(2) / prop.turning_radius.array() ).cast<double>() );
	
	if(I.size() != 0 ) {
		MatrixXd temp = vmag(I.reshaped(),0).array().pow(2).array() / prop.turning_radius(I.reshaped(),0).array();
		a_norm(I.reshaped(), all) = a_norm(I.reshaped(),all).array() / (a_norm_mag(I.reshaped(),0) * MatrixXd::Ones(1,3)).array() * (temp*MatrixXd::Ones(1,3)).array();
	}

	MatrixXd rates = append_right( states(all,seq(3,5)) , ( a_tan(all,seq(0,2)) + a_norm(all,seq(0,2)) ) );

	return rates;
}

// Apply velocity constraints
void sph_sim::constrain_vel() {
	// Max velocity constrain
	MatrixXd V = states(all,seq(3,5)).array().pow(2).rowwise().sum().sqrt();
	MatrixXi I = find( (V.array() > prop.vmax.array() ).cast<double>() );
	if(I.size() != 0 ) {
		states(I.reshaped(), 3) = states(I.reshaped(), 3).array() / V(I.reshaped(),0).array() * prop.vmax(I.reshaped(),0).array();
		states(I.reshaped(), 4) = states(I.reshaped(), 4).array() / V(I.reshaped(),0).array() * prop.vmax(I.reshaped(),0).array();
		states(I.reshaped(), 5) = states(I.reshaped(), 5).array() / V(I.reshaped(),0).array() * prop.vmax(I.reshaped(),0).array();
	}

	// Min velocity constrain
	V = states(all,seq(3,5)).array().pow(2).rowwise().sum().sqrt();
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
// Return the current time in the SPH simulation
double sph_sim::get_time() {
	return t;
}
// Return the initial time for the SPH simulation
double sph_sim::get_initial_time() {
	return t0;
}
// Return the time step to be used for the SPH simulation
double sph_sim::get_dt() {
	return param.dt;
}

// Return a matrix containing the [x y z] positions and [u v w] velocities of all SPH particles.
// Each particle is stored in one row:
//			[ x0 y0 z0 u0 v0 w0 ]
// states = [ x1 y1 z1 u1 v1 w1 ]
//			[		 ...		]
MatrixXd sph_sim::get_states() {
	return states;
}

// Return the total number of particles in the simulation
int sph_sim::get_npart() {
	return npart;
}
// Return the number of vehicles in the simulation
int sph_sim::get_nveh() {
	return nveh;
}
// Return the number of obstacles in the simulation
int sph_sim::get_nobs() {
	return nobs;
}
// Return the number of reduced density (attractor) particles in the simulation
int sph_sim::get_nrd() {
	return nrd;
}

// Return a column vector containing x/y/z positions or u/v/w velocities of all the SPH particles
MatrixXd sph_sim::get_x() {
	return states.col(0);
}
MatrixXd sph_sim::get_y() {
	return states.col(1);
}
MatrixXd sph_sim::get_z() {
	return states.col(2);
}
MatrixXd sph_sim::get_u() {
	return states.col(3);
}
MatrixXd sph_sim::get_v() {
	return states.col(4);
}
MatrixXd sph_sim::get_w() {
	return states.col(5);
}

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
prop_struct sph_sim::get_prop() {
	return prop;
}

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
param_struct sph_sim::get_param() {
	return param;
}

// Return a group_conf_struct containing the group configuration
group_conf_struct sph_sim::get_group_conf() {
	return group_conf;
}

// TODO
// INCOMPLETE, MODIFIED
// contemplate the type of vseq. int or double? eh...