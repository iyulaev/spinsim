/* 	llg1_6.c - part of Spinsim_GSL Project
	This file only contains the function llg1_6()
	It is used by spinsim.c as the differential equation to solve for the magnetization
	of a 3-layer structure, having one pinned and two free layers 
	Written and Engineered 2010-01-16 S. Moyerman
	University of California San Diego, ECE Department, Fullerton Lab
*/

//#define FMR_FIELD

//Turns on debugging output; won't run normally but will print out the fields involved
//#define LLG6_DEBUG_ON

//Removed mutual demag field from the simulation. Can be useful for simplifying simulation and isolating effects
//#define REMOVE_MUTUAL_DEMAG

//Swaps Hz[k] value and J_1,2 exchange parameter. Instead of sweeping field we are sweeping the J_1,2 exchange parameters
//#define LLG6_SWAP_JEXCH_AND_HK

//When the above define is enabled, this sets the DC field applied, so that a field can still be used. ONLY WORKS WHEN LLG6_SWAP_JEXCH_AND_HK IS ENABLED!!!
//#define HZ_K_VAL -1000000

//The below define, when un-commented, swaps the current density values with the J_1,2 exchange coefficient
//and sets the current density to 0
//#define LLG6_SWAP_JEXCH_AND_J

/* llg1_6() is a function that is the "differential equation" for the ODE solver, when solving for a system with three layers
It takes parameters as defined by GSL odeiv system (http://www.gnu.org/software/gsl/manual/html_node/Defining-the-ODE-System.html)
t = the time position to use
y[] = the values of y at which our equation is being computed±
dydt[] = the array in which to store the output
params = a pointer to a memory space of 23 pointers. See code below for how they are unpacked.
The function will return the values of dy/dt for each of the three dimensions STORED INTO dydt[]
The actual return value is GSL_SUCCESS for a valid execution */
int llg1_6(double t, const double y[], double dydt[], void * params) {
//	First, extract all values from params parameter
//	Ivan: We make a poor non-portable assumption that pointers to all types have same length
//	But it's probably good enough for most architectures

	double ** param_list = (double **) params;
	double * init_mag = param_list[0];
	double * Hz = param_list[1];
	double * j = param_list[2];
	double Beta_div_P_and_j = *param_list[3];
	double * Alpha = param_list[4];
	double * Alpha_bar = param_list[5];
	double* H_K1 = param_list[6];
	double* H_K2 = param_list[7];
	double thetaoffset = *param_list[8];
	gsl_matrix** ND = (gsl_matrix**) param_list[9];
	gsl_matrix** NR = (gsl_matrix**) param_list[10];
	gsl_matrix** NR_down = (gsl_matrix**) param_list[11];
	double* Ms = param_list[12];
	double * Gamma_bar = param_list[13];
	int k = *((int*)param_list[14]);
	int i = *((int*)param_list[15]);
	int j_timedep_funs = *((int*)param_list[16]);
	int hz_timedep_funs = *((int*)param_list[17]);
	double temperature = *param_list[18];
	double dx = *param_list[19];
	double dy = *param_list[20];
	double* dz = param_list[21];
	double* layer_position_z = param_list[22];
	double* af_coupling_js = param_list[23];
	double* pol = param_list[24];
	double t_relax = *param_list[25];
	int sim_field_in_plane = *((int *)param_list[26]);
	
	int count;
	double J_i = j[i];
	double Hz_k = Hz[k];
	
	#ifdef LLG6_SWAP_JEXCH_AND_HK
	double* temp_af_coupling_js = (double*) malloc(6 * sizeof(double));
	for(count = 0; count < 6; count++) temp_af_coupling_js[count] = af_coupling_js[count];
	temp_af_coupling_js[0] = Hz[k];
	
	Hz_k = HZ_K_VAL; //We'll use this in place of Hz[k] for the effective field applied externally
	af_coupling_js = temp_af_coupling_js; //Use temp_af_coupling_js instead of the one provided (so we don't kill it for everyone else)
	#endif
	
	#ifdef LLG6_SWAP_JEXCH_AND_J
	//Use temp_af_coupling_js  for coupling value storage instead of the one provided 
	//(so we don't kill it for everyone else)
	double* temp_af_coupling_js = (double*) malloc(6 * sizeof(double));
	for(count = 0; count < 6; count++) temp_af_coupling_js[count] = af_coupling_js[count];
	temp_af_coupling_js[0] = j[i];
	
	J_i = 0; //We'll use this in place of j[i] for the applied current density
	af_coupling_js = temp_af_coupling_js; 
	#endif
	
//	This is the current density we'll use for this iterations - pre-calculate it so we don't have to recalculate each time we use it
	double this_iterations_j = TimeDepFun(CALC_CURRENT, t, J_i, j_timedep_funs, t_relax);

//	Now we'll define and calculate variable's we'll actually use
//	Calculate spin transfer torque coefficient from layer 1 to layer 0
	double Beta_div_P_and_j_10 = (GAMMA*HBAR)/(2*(dz[0]*1e-9)*MU0*Ms[0]*Q_E);
	double Beta = Beta_div_P_and_j_10*polarisationpp(y[0],y[1],y[2],y[3],y[4],y[5], pol[0], pol[1])*this_iterations_j;
	double Beta_bar1 = Alpha[0]*Beta/( 1 + pow(Alpha[0],2) );
	double Beta_bar2 = Beta/(1+pow(Alpha[0],2));

//	 Spin transfer tq from layer 0 to layer 1
	double Beta_div_P_and_j_01 = (GAMMA*HBAR)/(2*(dz[1]*1e-9)*MU0*Ms[1]*Q_E);
	double Betap = -1 * Beta_div_P_and_j_01*polarisationpp(y[3],y[4],y[5],y[0],y[1],y[2], pol[1], pol[0])*this_iterations_j;
	double Beta_bar1p = Alpha[1]*Betap/( 1 + pow(Alpha[1],2) );
	double Beta_bar2p = Betap/(1+pow(Alpha[1],2));


//	Re-define p as p_mat, a matrix with the values as different columns
	gsl_matrix* p_mat = gsl_matrix_calloc(3,1);
	gsl_matrix_set(p_mat, 0, 0, y[3]);
	gsl_matrix_set(p_mat, 1, 0, y[4]);
	gsl_matrix_set(p_mat, 2, 0, y[5]);
	
//	Re-define y as y_mat, a matrix with the values as different columns
	gsl_matrix* y_mat = gsl_matrix_calloc(3,1);
	gsl_matrix_set(y_mat, 0, 0, y[0]);
	gsl_matrix_set(y_mat, 1, 0, y[1]);
	gsl_matrix_set(y_mat, 2, 0, y[2]);

//	  Free Layer
//	Define H, the applied field, as a gsl_matrix, remember that it's (rows, columns)
	gsl_matrix* H = gsl_matrix_calloc(3,1);
	#ifdef FMR_FIELD	
	if(sim_field_in_plane == 1) {
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax) - H[k]);
		gsl_matrix_set(H, 1, 0, 0.0);
		gsl_matrix_set(H, 0, 0, H[k]);
	}
	else {
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax) - H[k]);
		gsl_matrix_set(H, 1, 0, 0.0);
		gsl_matrix_set(H, 2, 0, H[k]);
	}
	#else
	gsl_matrix_set(H, 0, 0, 0.0);
	gsl_matrix_set(H, 1, 0, 0.0);
	gsl_matrix_set(H, 2, 0, 0.0);
	
	if(sim_field_in_plane == 1) {
		#ifdef LLG6_SWAP_JEXCH_AND_HK
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz_k, hz_timedep_funs, t_relax));
		#else
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax));
		#endif
	}
	else {
		#ifdef LLG6_SWAP_JEXCH_AND_HK
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz_k, hz_timedep_funs, t_relax));
		#else
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax));
		#endif
	}
	
	#endif
//	Anisotropic field
	gsl_matrix* Hani = gsl_matrix_calloc(3,1);
	gsl_matrix_set(Hani, 0, 0, H_K1[0]*y[0]*sin(thetaoffset));
	gsl_matrix_set(Hani, 1, 0, 0.0);
	gsl_matrix_set(Hani, 2, 0, H_K1[0]*y[2]*cos(thetaoffset)+H_K2[0]*y[2]*( pow(y[0],2) + pow(y[1],2) ) );
//	Self-demag field "Hani"
	gsl_matrix* HD = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND[0], y_mat, 0, HD);
	gsl_matrix_scale(HD, (-1.0 * Ms[0]));
//	Mutual demag field "HR" (from layer 2 to layer 1)

	gsl_matrix* HR = gsl_matrix_calloc(3,1);
#ifndef REMOVE_MUTUAL_DEMAG
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR[0], p_mat, 0, HR);
	gsl_matrix_scale(HR, (-1.0 * Ms[1]));
#endif
	
//	Anti-)Ferromagnetic Coupling Effective Field
	gsl_matrix* Hafc = gsl_matrix_calloc(3,1);
	
	//gsl_matrix_set(Hafc, 2, 0, af_coupling_js[0] * (-1.0) * cos_angle_bw_two_vectors(y[0], y[1], y[2], init_mag[3], init_mag[4], init_mag[5]));
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(Hafc, 0, 0, y[3] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 1, 0, y[4] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 2, 0, y[5] * af_coupling_js[0] / (1e6 * Ms[0] * MU0 * dz[0]*1e-9));
	#else
	gsl_matrix_set(Hafc, 0, 0, y[3] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 1, 0, y[4] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 2, 0, y[5] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	#endif
	
	//How to create a random thermal "effective field" vector
	//Generate the random components
	double ht_comps0[3], ht_comps_len0;
	gsl_matrix* Ht0 = gsl_matrix_alloc(3,1);
	
	//Don't bother calculating if temperature is 0!
	if(temperature > 0.0) {
		for(count = 0; count < 3; count++) ht_comps0[count] = (1.0 * rand())/(1.0 * RAND_MAX);
		//Normalize said random components
		ht_comps_len0 = sqrt(pow(ht_comps0[0], 2) + pow(ht_comps0[1], 2) + pow(ht_comps0[2], 2));
		for(count = 0; count < 3; count++) ht_comps0[count] = ht_comps0[count] / ht_comps_len0;
		//Calculate the magnitude based on Boltzmann's constant etc, and multiply into components
		//For the dimension of the sample, we use dx * dy * (layer position of top layer + thickness of top layer)
		ht_comps_len0 = (temperature * K_BOLTZMANN) / (Ms[0] * (dx * 1e-9 * dy * 1e-9 * (dz[1]+layer_position_z[1]) * 1e-9));
		for(count = 0; count < 3; count++) ht_comps0[count] = ht_comps0[count] * ht_comps_len0;
		//Put all of random components into a matrix
		gsl_matrix_set(Ht0, 0, 0, ht_comps0[0]);
		gsl_matrix_set(Ht0, 1, 0, ht_comps0[1]);
		gsl_matrix_set(Ht0, 2, 0, ht_comps0[2]);
	}

	
//	Effective field is the sum of the anisotropy field (Hani), the self (HD)
//	and mutual (HR) demag. fields, and the external field (H)
	gsl_matrix* Heff = gsl_matrix_calloc(3,1); //calloc sets initial values to 0
	gsl_matrix_add(Heff, H);
	gsl_matrix_add(Heff, Hani);
	gsl_matrix_add(Heff, HR);
	gsl_matrix_add(Heff, HD);
	gsl_matrix_add(Heff, Hafc);
	if(temperature > 0.0) {
		gsl_matrix_add(Heff, Ht0);
	}

//	  Pinned Layer
	gsl_matrix* HPani = gsl_matrix_calloc(3,1);
	gsl_matrix_set(HPani, 0, 0, 0.0);
	gsl_matrix_set(HPani, 1, 0, 0.0);
	gsl_matrix_set(HPani, 2, 0, H_K1[1]*y[5]);
	
	gsl_matrix* HPD = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND[1], p_mat, 0, HPD);
	gsl_matrix_scale(HPD, (-1.0 * Ms[1]));
	
	//Mutual demag tensor
	//From layer 1 to layer 2, use NR_down
	gsl_matrix* HPR = gsl_matrix_calloc(3,1);
#ifndef REMOVE_MUTUAL_DEMAG
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR_down[0], y_mat, 0, HPR);
	gsl_matrix_scale(HPR, (-1.0 * Ms[0]));
#endif
	
	//(Anti-)Ferromagnetic Coupling Effective Field
	gsl_matrix* HPafc = gsl_matrix_calloc(3,1);
	
	//gsl_matrix_set(Hafc, 2, 0, af_coupling_js[0] * (-1.0) * cos_angle_bw_two_vectors(y[0], y[1], y[2], init_mag[3], init_mag[4], init_mag[5]));
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(HPafc, 0, 0, y[0] * af_coupling_js[0] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 1, 0, y[1] * af_coupling_js[0] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 2, 0, y[2] * af_coupling_js[0] / (1e6 * Ms[1] * dz[1]*1e-9 * MU0));
	#else
	gsl_matrix_set(HPafc, 0, 0, y[0] * af_coupling_js[0] / (Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 1, 0, y[1] * af_coupling_js[0] / (Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 2, 0, y[2] * af_coupling_js[0] / (Ms[1] * dz[1]*1e-9));
	#endif
	
	//How to create a random thermal "effective field" vector
	//Generate the random components
	double ht_comps1[3], ht_comps_len1;
	gsl_matrix* Ht1 = gsl_matrix_alloc(3,1);
	//int count;
	
	//Don't bother calculating if temperature is 0!
	if(temperature > 0.0) {
		for(count = 0; count < 3; count++) ht_comps1[count] = (1.0 * rand())/(1.0 * RAND_MAX);
		//Normalize said random components
		ht_comps_len1 = sqrt(pow(ht_comps1[0], 2) + pow(ht_comps1[1], 2) + pow(ht_comps1[2], 2));
		for(count = 0; count < 3; count++) ht_comps1[count] = ht_comps1[count] / ht_comps_len1;
		//Calculate the magnitude based on Boltzmann's constant etc, and multiply into components
		//For the dimension of the sample, we use dx * dy * (layer position of top layer + thickness of top layer)
		ht_comps_len1 = (temperature * K_BOLTZMANN) / (Ms[1] * (dx * 1e-9 * dy * 1e-9 * (dz[1]+layer_position_z[1]) * 1e-9));
		for(count = 0; count < 3; count++) ht_comps1[count] = ht_comps1[count] * ht_comps_len1;
		//Put all of random components into a matrix
		gsl_matrix_set(Ht1, 0, 0, ht_comps1[0]);
		gsl_matrix_set(Ht1, 1, 0, ht_comps1[1]);
		gsl_matrix_set(Ht1, 2, 0, ht_comps1[2]);
	}
	
	gsl_matrix* HPeff = gsl_matrix_calloc(3,1); //calloc sets the initival values to 0											 
	gsl_matrix_set_zero(HPeff);
	gsl_matrix_add(HPeff, H);
	gsl_matrix_add(HPeff, HPani);
	gsl_matrix_add(HPeff, HPR);
	gsl_matrix_add(HPeff, HPD);
	gsl_matrix_add(HPeff, HPafc);
	if(temperature > 0.0) {
		gsl_matrix_add(HPeff, Ht1);
	}

//	In fact, one has to integrate three differential equations,
//	one for each component of m

//The original matlab code
/*dy(1) = -Gamma_bar*(y(2)*Heff(3)-y(3)*Heff(2))-Alpha_bar*((y(1)*y(1)-1)*Heff(1)+...
		 y(1)*y(2)*Heff(2)+y(1)*y(3)*Heff(3))+Beta_bar1*(y(2)*y(6) - y(3)*y(5))...
		   -Beta_bar2*((y(1)*y(1)-1)*y(4)+y(1)*y(2)*y(5)+y(1)*y(3)*y(6));

	 dy(2) = -Gamma_bar*(y(3)*Heff(1)-y(1)*Heff(3))-Alpha_bar*((y(2)*y(2)-1)*Heff(2)+...
		 y(1)*y(2)*Heff(1)+y(2)*y(3)*Heff(3))+Beta_bar1*(y(3)*y(4)-y(1)*y(6))...
		   -Beta_bar2*((y(2)*y(2)-1)*y(5)+y(1)*y(2)*y(4)+y(2)*y(3)*y(6));

	 dy(3) = -Gamma_bar*(y(1)*Heff(2)-y(2)*Heff(1))-Alpha_bar*((y(3)*y(3)-1)*Heff(3)...
		 +y(1)*y(3)*Heff(1)+y(2)*y(3)*Heff(2))+Beta_bar1*(y(1)*y(5)-y(2)*y(4))....
		   -Beta_bar2*((y(3)*y(3)-1)*y(6)+y(1)*y(3)*y(4)+y(2)*y(3)*y(5));
	
	 dy(4) = -Gamma_bar*(y(5)*HPeff(3)-y(6)*HPeff(2))-Alpha_bar*((y(4)*y(4)-1)*HPeff(1)+...
		 y(4)*y(5)*HPeff(2)+y(4)*y(6)*HPeff(3))+Beta_bar1p*(y(5)*y(3) - y(6)*y(2))...
		   -Beta_bar2p*((y(4)*y(4)-1)*y(1)+y(4)*y(5)*y(2)+y(4)*y(6)*y(3));

	 dy(5) = -Gamma_bar*(y(6)*HPeff(1)-y(4)*HPeff(3))-Alpha_bar*((y(5)*y(5)-1)*HPeff(2)+...
		 y(4)*y(5)*HPeff(1)+y(5)*y(6)*HPeff(3))+Beta_bar1p*(y(6)*y(1)-y(4)*y(3))...
		   -Beta_bar2p*((y(5)*y(5)-1)*y(2)+y(4)*y(5)*y(1)+y(5)*y(6)*y(3));

	 dy(6) = -Gamma_bar*(y(4)*HPeff(2)-y(5)*HPeff(1))-Alpha_bar*((y(6)*y(6)-1)*HPeff(3)...
		 +y(4)*y(6)*HPeff(1)+y(5)*y(6)*HPeff(2))+Beta_bar1p*(y(4)*y(2)-y(5)*y(1))....
		   -Beta_bar2p*((y(6)*y(6)-1)*y(3)+y(4)*y(6)*y(1)+y(5)*y(6)*y(2));*/

	// Free Layer Dynamics
	dydt[0] = -1.0 * Gamma_bar[0] * (y[1]*gsl_matrix_get(Heff, 2, 0)-y[2]*gsl_matrix_get(Heff, 1, 0)) - Alpha_bar[0]*((y[0]*y[0]-1)*gsl_matrix_get(Heff, 0, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 1, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 2, 0)) + Beta_bar1 * (y[1]*y[5] - y[2]*y[4]) - Beta_bar2 * ((y[0]*y[0]-1)*y[3]+y[0]*y[1]*y[4]+y[0]*y[2]*y[5]);
	dydt[1] = -1.0 * Gamma_bar[0]*(y[2]*gsl_matrix_get(Heff, 0, 0)-y[0]*gsl_matrix_get(Heff, 2, 0))-Alpha_bar[0]*((y[1]*y[1]-1)*gsl_matrix_get(Heff, 1, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 2, 0))+Beta_bar1*(y[2]*y[3]-y[0]*y[5]) - Beta_bar2*((y[1]*y[1]-1)*y[4]+y[0]*y[1]*y[3]+y[1]*y[2]*y[5]);
	dydt[2] = -1.0 * Gamma_bar[0]*(y[0]*gsl_matrix_get(Heff, 1, 0)-y[1]*gsl_matrix_get(Heff, 0, 0))-Alpha_bar[0]*((y[2]*y[2]-1)*gsl_matrix_get(Heff, 2, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 1, 0))+Beta_bar1*(y[0]*y[4]-y[1]*y[3])-Beta_bar2*((y[2]*y[2]-1)*y[5]+y[0]*y[2]*y[3]+y[1]*y[2]*y[4]);
	
	// Ref Layer Dynamics
	dydt[3] = -1.0 * Gamma_bar[1]*(y[4]*gsl_matrix_get(HPeff, 2, 0)-y[5]*gsl_matrix_get(HPeff, 1, 0)) - Alpha_bar[1]*((y[3]*y[3]-1)*gsl_matrix_get(HPeff, 0, 0)+y[3]*y[4]*gsl_matrix_get(HPeff, 1, 0)+y[3]*y[5]*gsl_matrix_get(HPeff, 2, 0))+Beta_bar1p*(y[4]*y[2]-y[5]*y[1])-Beta_bar2p*((y[3]*y[3]-1)*y[0]+y[3]*y[4]*y[1]+y[3]*y[5]*y[2]);
	dydt[4] = -1.0*Gamma_bar[1]*(y[5]*gsl_matrix_get(HPeff, 0, 0)-y[3]*gsl_matrix_get(HPeff, 2, 0))-Alpha_bar[1]*((y[4]*y[4]-1)*gsl_matrix_get(HPeff, 1, 0)+y[3]*y[4]*gsl_matrix_get(HPeff, 0, 0)+y[4]*y[5]*gsl_matrix_get(HPeff, 2, 0))+Beta_bar1p*(y[5]*y[0]-y[3]*y[2])-Beta_bar2p*((y[4]*y[4]-1)*y[1]+y[3]*y[4]*y[0]+y[4]*y[5]*y[2]);
	dydt[5] = -1.0*Gamma_bar[1]*(y[3]*gsl_matrix_get(HPeff, 1, 0)-y[4]*gsl_matrix_get(HPeff, 0, 0))-Alpha_bar[1]*((y[5]*y[5]-1)*gsl_matrix_get(HPeff, 2, 0)+y[3]*y[5]*gsl_matrix_get(HPeff, 0, 0)+y[4]*y[5]*gsl_matrix_get(HPeff, 1, 0))+Beta_bar1p*(y[3]*y[1]-y[4]*y[0])-Beta_bar2p*((y[5]*y[5]-1)*y[2]+y[3]*y[5]*y[0]+y[4]*y[5]*y[1]);
	
	//debug statements
#ifdef LLG6_DEBUG_ON
		//fprintf(stderr,"H_K1[0] = %e, H_K2[0] = %e\n\n", H_K1[0], H_K2[0]);
		fprintf(stderr,"Heff = %e %e %e\n\n", gsl_matrix_get(Heff, 0, 0), gsl_matrix_get(Heff, 1, 0), gsl_matrix_get(Heff, 2, 0));
		//ffprintf(stderr,stderr, "Got here 1\n");
		fprintf(stderr,"Hani = %e %e %e\n\n", gsl_matrix_get(Hani, 0, 0), gsl_matrix_get(Hani, 1, 0), gsl_matrix_get(Hani, 2, 0));
		fprintf(stderr,"H = %e %e %e\n\n", gsl_matrix_get(H, 0, 0), gsl_matrix_get(H, 1, 0), gsl_matrix_get(H, 2, 0));
		fprintf(stderr,"HR = %e %e %e\n\n", gsl_matrix_get(HR, 0, 0), gsl_matrix_get(HR, 1, 0), gsl_matrix_get(HR, 2, 0));
		fprintf(stderr,"HD = %e %e %e\n\n", gsl_matrix_get(HD, 0, 0), gsl_matrix_get(HD, 1, 0), gsl_matrix_get(HD, 2, 0));
		fprintf(stderr,"Hafc = %e %e %e\n\n", gsl_matrix_get(Hafc, 0, 0), gsl_matrix_get(Hafc, 1, 0), gsl_matrix_get(Hafc, 2, 0));
		
		fprintf(stderr,"HPeff = %e %e %e\n\n", gsl_matrix_get(HPeff, 0, 0), gsl_matrix_get(HPeff, 1, 0), gsl_matrix_get(HPeff, 2, 0));
		//ffprintf(stderr,stderr, "Got here 1\n");
		fprintf(stderr,"HPani = %e %e %e\n\n", gsl_matrix_get(HPani, 0, 0), gsl_matrix_get(HPani, 1, 0), gsl_matrix_get(HPani, 2, 0));
		fprintf(stderr,"HPR = %e %e %e\n\n", gsl_matrix_get(HPR, 0, 0), gsl_matrix_get(HPR, 1, 0), gsl_matrix_get(HPR, 2, 0));
		fprintf(stderr,"HPD = %e %e %e\n\n", gsl_matrix_get(HPD, 0, 0), gsl_matrix_get(HPD, 1, 0), gsl_matrix_get(HPD, 2, 0));
		fprintf(stderr,"HPafc = %e %e %e\n\n", gsl_matrix_get(HPafc, 0, 0), gsl_matrix_get(HPafc, 1, 0), gsl_matrix_get(HPafc, 2, 0));
		
		//fprintf(stderr,"Gamma_bar = %e, Alpha_bar = %e, y[] = (%e, %e, %e), init_mag[] = (%e, %e, %e)\n\n", Gamma_bar, Alpha_bar, y[0], y[1], y[2], init_mag[3], init_mag[4], init_mag[5]);
		fprintf(stderr,"Alpha[0] = %e, Alpha[1] = %e\n", Alpha[0], Alpha[1]);
		fprintf(stderr,"Alpha_bar[0] = %e, Alpha_bar[1] = %e\n", Alpha_bar[0], Alpha_bar[1]);
		/*fprintf(stderr,"Beta_bar1 = %e, Beta_bar2 = %e\n", Beta_bar1, Beta_bar2);
		fprintf(stderr,"Beta_div_P_and_j = %e\n", Beta_div_P_and_j);
		fprintf(stderr,"dYdT = (%e, %e, %e)\t(%e, %e, %e)\n", dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5]);*/
		printf("100%\n");
		exit(0);
#endif
	
	//free things we've used, mostly all of the matrices
	
	//free free layer
	gsl_matrix_free(Heff);
	gsl_matrix_free(HR);	
	gsl_matrix_free(Hani);
	gsl_matrix_free(HD);
	gsl_matrix_free(H);
	
	//free pinned layer
	gsl_matrix_free(HPani);
	gsl_matrix_free(HPeff);
	gsl_matrix_free(HPR);
	gsl_matrix_free(HPD);
	
	gsl_matrix_free(p_mat);
	gsl_matrix_free(y_mat);
	
	gsl_matrix_free(Hafc);
	gsl_matrix_free(HPafc);
	
	gsl_matrix_free(Ht0);
	gsl_matrix_free(Ht1);
	
	#ifdef LLG6_SWAP_JEXCH_AND_HK
	free(temp_af_coupling_js);
	#endif
	
	#ifdef LLG6_SWAP_JEXCH_AND_J
	free(temp_af_coupling_js);
	#endif
	
	return(GSL_SUCCESS);
}

