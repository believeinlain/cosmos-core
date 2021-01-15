#include "sph_sim.h"

// Evaluates the smoothing kernel function for the SPH equations
/*matrix kernel(matrix r, double h, int type) {
	matrix s = r/h;
	
	matrix W(s.size());

	// Cubic spline kernel, used for vehicle-reduced density (type 1) particle interactions.
	//W = W + (type==1).*( ( 1 - 3/2*s.^2        + 3/4*s.^3 )          .*( s<1 ) + ( 1/4*(2-s).^3 )           .*( (s >= 1)   .*(s <= 2) ) )   ./(pi*h.^3);
	W = W + (type==1)*( ( 1.0-3.0/2.0*pow(s,2.0) + 3.0/4.0*pow(s,3.0) ) * (s<1) + ( 1.0/4.0*pow((2.0-s),3.0) ) * ( (s >= 1.0) * (s <= 2.0) ) ) / M_PI*pow(h,3.0);

	// Quadratic kernel, used for all other (type 2) interactions
	//W = W + ( type==2 ).*15./(16*pi*h.^3)     .*(s.^2/4-s+1)         .*(s<2);
	W = W + (type==1)*15.0/(16.0*M_PI*pow(h,3.0))*(pow(s,2.0),4.0-s+1.0)*(s<2.0);
	return W;
}*/
// Evaluates the smoothing kernel function for the SPH equations
double kernel(double r, double h, int type) {
	double s = r/h;
	
	double W = 0;


	// Cubic spline kernel, used for vehicle-reduced density (type 1) particle interactions.
	//W = W + (type==1).*( ( 1 - 3/2*s.^2        + 3/4*s.^3 )          .*( s<1 ) + ( 1/4*(2-s).^3 )           .*( (s >= 1)   .*(s <= 2) ) )   ./(pi*h.^3);
	W = W + (type==1)*( ( 1.0-3.0/2.0*pow(s,2.0) + 3.0/4.0*pow(s,3.0) ) * (s<1) + ( 1.0/4.0*pow((2.0-s),3.0) ) * ( (s >= 1.0) * (s <= 2.0) ) ) / M_PI*pow(h,3.0);

	// Quadratic kernel, used for all other (type 2) interactions
	//W = W + ( type==2 ).*15./(16*pi*h.^3)     .*(s.^2/4-s+1)         .*(s<2);
	W = W + (type==1)*15.0/(16.0*M_PI*pow(h,3.0))*(pow(s,2.0),4.0-s+1.0)*(s<2.0);
	return W;
}

// Evaluates the derivatie of the smoothing kernel function for the SPH equations.
// Note that this returns a scalar value which is the magnitude of the gradient W.
// The direction (if needed) must be computed separately.
double kernel_grad(double r, double h, int type) {
	double s = r/h;

	double dWdr = 0;

	// Cubic spline kernel, used for vehicle-reduced density (type 1) particle interactions.
	//dWdr = dWdr + ( type==1 ).*( ( -3*s + 9/4*s.^2 )     .*(s<1) + ( -3/4*(2-s).^2 )           .*( (s >= 1).*(s <= 2) ) )./(pi*h.^4);
	dWdr = dWdr + (type==1)*( ( -3.0*s + 9.0/4.0*pow(s,2.0))*(s<1) + ( -3.0/4.0*pow((2.0-s),2.0) )*( (s>=1.0)*(s<=2.0)  ) )/(M_PI*pow(h,4.0));

	// Quadratic kernel, used for all other (type 2) interactions
	//dWdr = dWdr + ( type==2 ).*15./(16*pi*h.^4)     .*(s/2-1)   .*(s<2);
	dWdr = dWdr + (type==2)*15.0/(16.0*M_PI*pow(h,4.0))*(s/2.0-1.0)*(s<2.0);

	return dWdr;
}

sph_sim::sph_sim(param_struct param, group_conf_struct group_conf, double t0 = 0) {
	sph_sim::param = param;
	sph_sim::group_conf = group_conf;
	t0 = t0;

	// Setup SPH properties
	init_prop();
	compute_hij();
	kernel_type();

	// Initialize positions and velocities
	init_states();
	
}

void sph_sim::init_prop() {
	int N = 0;

	// Initialize vehicles
	for(int i = 0; i < this->group_conf.num_veh.size(); ++i) {
		for(int j = 0; j < this->group_conf.num_veh[i]; ++j) {
			// Motion constraints
			this->prop.vmin(N,0) = this->group_conf.veh_limits.vmin(i);
			this->prop.vmax(N,0) = this->group_conf.veh_limits.vmax(i);
			this->prop.turning_radius(N,0) = this->group_conf.veh_limits.turning_radius(i);
			this->prop.amax(N,0) = pow(this->prop.vmax(N),2.0) / this->prop.turning_radius(N);

			// Smoothing width
			this->prop.h(N,0) = this->group_conf.veh_h(i);

			// Mass
			this->prop.m(N,0) = rho0 / kernel(0, this->prop.h(N), 2);

			// Kernel values at 0 and h
			double KER0 = kernel(0, this->prop.h(N), 2);
			double KERh = kernel(this->prop.h(N), this->prop.h(N), 2);

			// Kernel gradient at h
			double KERG = kernel_grad(this->prop.h(N), this->prop.h(N), 2);

			// Pressure force ~ K*Fp
			double Fp = rho0 * KERh / KER0;
			// Viscous force ~ mu*Fmu
			double Fmu = 2.0 * this->prop.vmax(N) * KER0 * KERG / ( rho0 * this->prop.h(N) * pow((KER0+KERh),2.0) );

			// Force coeficients found by solving:
			//	amax = vmax^2/turning_radius = K*Fp + mu*Fmu
			//	Re = (K*Fp)/(mu*Fmu)
			//
			// This enforces the desired Reynolds number at r_ij=h
			// and limits the acceleration magnitude to be on the
			// order of amax.
			this->prop.mu(N,0) = this->prop.vmax(N) / this->prop.turning_radius(N) / Fmu / (1.0 + this->param.Re);
			this->prop.K(N,0) = this->param.accel.veh * this->param.Re * this->prop.mu(N) * Fmu / Fp;

			// Group number and particle type
			this->prop.group(N,0) = i;
			this->prop.particle_type(N,0) = particle_type_enum::veh;	// NOTE: may produce oob error?
			N++;
		}
	}

	this->nveh = N;

	// Obstacles
	for(int i = 0; i < this->group_conf.num_obs; ++i) {
		// Motion constraints
		this->prop.vmin(N,0) = 0;
		this->prop.vmax(N,0) = 0;
		this->prop.turning_radius(N,0) = 0;
		this->prop.amax(N,0) = 0;

		// Smoothing width
		this->prop.h(N,0) = this->group_conf.obs_h(i);

		// Mass
		this->prop.m(N,0) = 2.0 * rho0 / kernel(0, this->prop.h(N), 2);

		// Kernel values at 0 and h
		double KER0 = kernel(0, this->prop.h(N), 2);
		double KERh = kernel(this->prop.h(N), this->prop.h(N), 2);

		// Kernel gradient at h
		double KERG = kernel_grad(this->prop.h(N), this->prop.h(N), 2);

		// Force coeficients found by solving:
		this->prop.mu(N,0) = 0;
		this->prop.K(N,0) = this->param.accel.obs * this->prop.amax.max() * KER0 / (rho0 * KERh);

		// Group number and particle type
		this->prop.group(N,0) = i;
		this->prop.particle_type(N,0) = particle_type_enum::obs;
		N++;
	}

	this->nobs = this->group_conf.num_obs;

	// Reduced density particles
	for(int i = 0; i < this->group_conf.num_rd; ++i) {
		// Motion constraints
		this->prop.vmin(N,0) = 0;
		this->prop.vmax(N,0) = 0;
		this->prop.turning_radius(N,0) = 0;
		this->prop.amax(N,0) = 0;

		// Smoothing width
		this->prop.h(N,0) = this->group_conf.rd_h(i);

		// Mass
		this->prop.m(N,0) = rho0 / kernel(0, this->prop.h(N), 1) * 1e-8;

		// Force coeficients found by solving:
		// No viscosity for attractors
		this->prop.mu(N,0) = 0;
		matrix I = find(this->prop.group(range(0,this->nveh)) == this->group_conf.rd_group(i));
		this->prop.K(N,0) = -1.0 * this->param.accel.rd * this->prop.amax(I).max()
							* kernel(0, this->prop.h(N),1) / kernel_grad(this->prop.h(N),this->prop.h(N), 1);

		// Group number and particle type
		this->prop.group(N,0) = i;
		this->prop.particle_type(N,0) = particle_type_enum::obs;
		N++;
	}
}

// Compute h_ij matrix
void sph_sim::compute_hij() {
	matrix hi = this->prop.h(range(0,this->nveh)) * ones(1,this->nveh);
	matrix hj = hi.T();

	// Vehicles: hij=max(hi,hj)
	this->prop.hij = max(hi,hj);

	// Obstacles: hij=h_obs
	int I1 = this->nveh+1; // NOTE: +1 for 1 based indexing? Maybe not necessary. Consider I2 as well
	int I2 = this->nveh + this->group_conf.num_obs;
	this->prop.hij.append_down( this->prop.h(range(I1,I2)) * ones(1,this->prop.hij.dim().second) );
	this->prop.hij.append_right( ones(this->prop.hij.dim().first,1) * this->prop.h(range(I1,I2)).T() );

	// Reduced density particles: hij=h_rd
	I1 = this->nveh + this->group_conf.num_obs + 1;
	I2 = this->nveh + this->group_conf.num_obs + this->group_conf.num_rd;
	this->prop.hij.append_down( this->prop.h(range(I1,I2)) * ones(1,this->prop.hij.dim().second) );
	this->prop.hij.append_right( ones(this->prop.hij.dim().first,1) * this->prop.h(range(I1,I2)).T() );
}

// Create a matrix kernel_type that tells which kernel to use.
// 1 is for vehicle-reduced density particle interactions,
// 2 is for all others
void sph_sim::kernel_type() {
	int N = this->prop.m.dim().first;
	matrix ki = this->prop.particle_type * ones(1,N);
	matrix kj = ki.T();

	this->prop.kernel_type = 2*ones(N);

	matrix I = find(mlogical_and( (ki == particle_type_enum::veh * ones(N)) , (kj == particle_type_enum::rd * ones(N)) ));
	this->prop.kernel_type.bulk_assign(I,1);
	I = find(mlogical_and( (kj == particle_type_enum::veh * ones(N)) , (ki == particle_type_enum::rd * ones(N)) ));
	this->prop.kernel_type.bulk_assign(I,1);
}

// Set the initial SPH states (positions and velocities) for all particles
void sph_sim::init_states() {
	if(this->param.ndim == 2) {
		// 2d initialization
		this->init2d();
	} else {
		// 3d initialization
		this->init3d();
	}

	// Obstacles
	for(size_t i = 0; i < this->group_conf.num_obs; ++i) {
		this->states.assign_row(this->states.dim().first, { this->group_conf.obs_init.x[i],
															this->group_conf.obs_init.y[i],
															this->group_conf.obs_init.z[i],
															0, 0, 0 });	// u v w
	}

	// Reduced density particles
	this->states.append_down(matrix(this->group_conf.num_rd,6));

}


// NOTES
// when indexing by stuff like num_obs, consider that c++ is zero-based. Also, my range function is inclusive. May cause problems.