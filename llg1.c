/* 	llg1.c - part of Spinsim_GSL Project
	This file only contains the function llg1()
	It is used by spinsim.c as the differential equation to solve for the magnetization
	of a 2-layer structure 
	Written and engineered S. Moyerman
	University of California San Diego, ECE Department, Fullerton Lab
*/

//#define LLG1_DEBUG_ON
//#define FMR_FIELD

//Uncomment the below to remove self & mutual dipole
//#define LLG1_NO_SELFDIPOLE
//#define LLG1_NO_MUTDIPOLE

/* llg1() is a function that is the "differential equation" for the ODE solver
It takes parameters as defined by GSL odeiv system (http://www.gnu.org/software/gsl/manual/html_node/Defining-the-ODE-System.html)
t = the time position to use
y[] = the values of y at which our equation is being computed
dydt[] = the array in which to store the output
params = a pointer to a memory space of 23 pointers. See code below for how they are unpacked.
The function will return the values of dy/dt for each of the three dimensions STORED INTO dydt[]
The actual return value is GSL_SUCCESS for a valid execution */
int llg1(double t, const double y[], double dydt[], void * params) {
//	First, extract all values from params parameter
//	Ivan: We make a poor non-portable assumption that pointers to all types have same length
//	But it's probably good enough for most architectures
//	For what all of these parameters actually mean, please see the main() function in spinsim.c
//	All of the parameter names should be the same, or very similar
	double ** param_list = (double **) params;
	double * init_mag = param_list[0];
	double * Hz = param_list[1];
	double * j = param_list[2];
	double Beta_div_P_and_j = *param_list[3];
	double Alpha = (param_list[4])[0];
	double Alpha_bar = (param_list[5])[0];
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
	
//	Now we'll define and calculate variable's we'll actually use
//	We assume that the pinned layer is layer 1, and free layer is layer 0
	double Beta = Beta_div_P_and_j*polarisationpp(y[0],y[1],y[2],init_mag[3], init_mag[4], init_mag[5], pol[0], pol[1])*TimeDepFun(CALC_CURRENT, t, j[i], j_timedep_funs, t_relax);
	double Beta_bar1 = Alpha*Beta/( 1 + pow(Alpha,2) );
	double Beta_bar2 = Beta/(1+pow(Alpha,2));
	
	//How to create a random thermal "effective field" vector
	//Generate the random components
	double ht_comps[3], ht_comps_len;
	gsl_matrix* Ht = gsl_matrix_alloc(3,1);
	int count;
	
	//Don't bother calculating if temperature is 0!
	if(temperature > 0.0) {
		for(count = 0; count < 3; count++) ht_comps[count] = (1.0 * rand())/(1.0 * RAND_MAX);
		//Normalize said random components
		ht_comps_len = sqrt(pow(ht_comps[0], 2) + pow(ht_comps[1], 2) + pow(ht_comps[2], 2));
		for(count = 0; count < 3; count++) ht_comps[count] = ht_comps[count] / ht_comps_len;
		//Calculate the magnitude based on Boltzmann's constant etc, and multiply into components
		//For the dimension of the sample, we use dx * dy * (layer position of top layer + thickness of top layer)
		ht_comps_len = (temperature * K_BOLTZMANN) / (Ms[0] * (dx * 1e-9 * dy * 1e-9 * (dz[1]+layer_position_z[1]) * 1e-9));
		for(count = 0; count < 3; count++) ht_comps[count] = ht_comps[count] * ht_comps_len;
		//Put all of random components into a matrix
		gsl_matrix_set(Ht, 0, 0, ht_comps[0]);
		gsl_matrix_set(Ht, 1, 0, ht_comps[1]);
		gsl_matrix_set(Ht, 2, 0, ht_comps[2]);
	}

	
//	Re-define y as y_mat, a matrix with the values as different columns
	gsl_matrix* y_mat = gsl_matrix_alloc(3,1);
	gsl_matrix_set(y_mat, 0, 0, y[0]);
	gsl_matrix_set(y_mat, 1, 0, y[1]);
	gsl_matrix_set(y_mat, 2, 0, y[2]);
//	Re-define p as p_mat, a matrix with the values as different columns
	gsl_matrix* p_mat = gsl_matrix_alloc(3,1);
	gsl_matrix_set(p_mat, 0, 0, init_mag[3]);
	gsl_matrix_set(p_mat, 1, 0, init_mag[4]);
	gsl_matrix_set(p_mat, 2, 0, init_mag[5]);
//	Define H, the applied field, as a gsl_matrix, remember that it's (rows, columns)
	gsl_matrix* H = gsl_matrix_calloc(3,1);
	#ifdef FMR_FIELD
	if(sim_field_in_plane == 1) {
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax) - Hz[k]);
		gsl_matrix_set(H, 1, 0, 0.0);
		gsl_matrix_set(H, 0, 0, Hz[k]);
	}
	else {
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax) - Hz[k]);
		gsl_matrix_set(H, 1, 0, 0.0);
		gsl_matrix_set(H, 2, 0, Hz[k]);
	}
	#else
	if(sim_field_in_plane == 1) {
		gsl_matrix_set(H, 2, 0, 0.0);
		gsl_matrix_set(H, 1, 0, 0.0);
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax));
	}
	else {
		gsl_matrix_set(H, 0, 0, 0.0);
		gsl_matrix_set(H, 1, 0, 0.0);
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax));
	}
	#endif
//	Anisotropic field
	gsl_matrix* Hani = gsl_matrix_alloc(3,1);
	gsl_matrix_set(Hani, 0, 0, H_K1[0]*y[0]*sin(thetaoffset));
	gsl_matrix_set(Hani, 1, 0, 0.0);
	gsl_matrix_set(Hani, 2, 0, H_K1[0]*y[2]*cos(thetaoffset)+H_K2[0]*y[2]*( pow(y[0],2) + pow(y[1],2) ) );
//	Self-demag field "HD"
	gsl_matrix* HD = gsl_matrix_calloc(3,1);
	#ifdef LLG1_NO_SELFDIPOLE
	#else
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND[0], y_mat, 0, HD);
	gsl_matrix_scale(HD, (-1.0 * Ms[0]));
	#endif
//	Mutual demag field "HR", from layer 2 to layer 1
	gsl_matrix* HR = gsl_matrix_calloc(3,1);
	#ifdef LLG1_NO_MUTDIPOLE
	#else
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR[0], p_mat, 0, HR);
	gsl_matrix_scale(HR, (-1.0 * Ms[1]));
	#endif
//	(Anti-)Ferromagnetic Coupling Effective Field
	gsl_matrix* Hafc = gsl_matrix_alloc(3,1);
	
	//gsl_matrix_set(Hafc, 2, 0, af_coupling_js[0] * (-1.0) * cos_angle_bw_two_vectors(y[0], y[1], y[2], init_mag[3], init_mag[4], init_mag[5]));
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(Hafc, 0, 0, init_mag[3] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 1, 0, init_mag[4] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 2, 0, init_mag[5] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	#else
	gsl_matrix_set(Hafc, 0, 0, init_mag[3] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 1, 0, init_mag[4] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 2, 0, init_mag[5] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	#endif
//	Effective field is the sum of the anisotropy field (Hani), the self (HD)
//	and mutual (HR) demag. fields, and the external field (H)
	gsl_matrix* Heff = gsl_matrix_calloc(3,1); //calloc sets initial values to 0
	gsl_matrix_add(Heff, H);
	gsl_matrix_add(Heff, Hani);
	gsl_matrix_add(Heff, HR);
	gsl_matrix_add(Heff, HD);
	gsl_matrix_add(Heff, Hafc);

	if(temperature > 0.0) gsl_matrix_add(Heff, Ht);

	
//	In fact, one has to integrate three differential equations,
//	one for each component of m
//	Ivan: I don't understand this next line, but it was commented out anyway
//	dm_x/dt=fx(m_x,m_y,m_z,t)

//	The following is thus the analytical expression for each differential component of the field

	dydt[0] = -1.0 * Gamma_bar[0] * (y[1]*gsl_matrix_get(Heff, 2, 0)-y[2]*gsl_matrix_get(Heff, 1, 0)) - Alpha_bar*((y[0]*y[0]-1)*gsl_matrix_get(Heff, 0, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 1, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 2, 0)) + Beta_bar1 * (y[1]*init_mag[5] - y[2]*init_mag[4]) - Beta_bar2 * ((y[0]*y[0]-1)*init_mag[3]+y[0]*y[1]*init_mag[4]+y[0]*y[2]*init_mag[5]);
	dydt[1] = -1.0 * Gamma_bar[0] *(y[2]*gsl_matrix_get(Heff, 0, 0)-y[0]*gsl_matrix_get(Heff, 2, 0))-Alpha_bar*((y[1]*y[1]-1)*gsl_matrix_get(Heff, 1, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 2, 0))+Beta_bar1*(y[2]*init_mag[3]-y[0]*init_mag[5]) - Beta_bar2*((y[1]*y[1]-1)*init_mag[4]+y[0]*y[1]*init_mag[3]+y[1]*y[2]*init_mag[5]);
	dydt[2] = -1.0 * Gamma_bar[0]*(y[0]*gsl_matrix_get(Heff, 1, 0)-y[1]*gsl_matrix_get(Heff, 0, 0))-Alpha_bar*((y[2]*y[2]-1)*gsl_matrix_get(Heff, 2, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 1, 0))+Beta_bar1*(y[0]*init_mag[4]-y[1]*init_mag[3])-Beta_bar2*((y[2]*y[2]-1)*init_mag[5]+y[0]*y[2]*init_mag[3]+y[1]*y[2]*init_mag[4]);
	
	//debug statements
#ifdef LLG1_DEBUG_ON
		fprintf(stderr,"H_K1[0] = %e, H_K2[0] = %e\n\n", H_K1[0], H_K2[0]);
		fprintf(stderr,"Heff = %e %e %e\n\n", gsl_matrix_get(Heff, 0, 0), gsl_matrix_get(Heff, 1, 0), gsl_matrix_get(Heff, 2, 0));
		//ffprintf(stderr,stderr, "Got here 1\n");
		fprintf(stderr,"Hani = %e %e %e\n\n", gsl_matrix_get(Hani, 0, 0), gsl_matrix_get(Hani, 1, 0), gsl_matrix_get(Hani, 2, 0));
		fprintf(stderr,"H = %e %e %e\n\n", gsl_matrix_get(H, 0, 0), gsl_matrix_get(H, 1, 0), gsl_matrix_get(H, 2, 0));
		fprintf(stderr,"HR = %e %e %e\n\n", gsl_matrix_get(HR, 0, 0), gsl_matrix_get(HR, 1, 0), gsl_matrix_get(HR, 2, 0));
		fprintf(stderr,"HD = %e %e %e\n\n", gsl_matrix_get(HD, 0, 0), gsl_matrix_get(HD, 1, 0), gsl_matrix_get(HD, 2, 0));
		fprintf(stderr,"NR = %e %e %e\n\n", gsl_matrix_get(NR[0], 0, 0), gsl_matrix_get(NR[0], 1, 1), gsl_matrix_get(NR[0], 2, 2));
		fprintf(stderr,"ND = %e %e %e\n\n", gsl_matrix_get(ND[0], 0, 0), gsl_matrix_get(ND[0], 1, 1), gsl_matrix_get(ND[0], 2, 2));
//		fprintf(stderr,"Gamma_bar = %e, Alpha_bar = %e, y[] = (%e, %e, %e), init_mag[] = (%e, %e, %e)\n\n", Gamma_bar[0], Alpha_bar[0], y[0], y[1], y[2], init_mag[3], init_mag[4], init_mag[5]);
//		fprintf(stderr,"Alpha = %e\n", Alpha[0]);
//		fprintf(stderr,"Beta_bar1 = %e, Beta_bar2 = %e\n", Beta_bar1, Beta_bar2);
//		fprintf(stderr,"Beta_div_P_and_j = %e\n", Beta_div_P_and_j);
		fprintf(stderr,"dYdT = (%e, %e, %e)\n", dydt[0], dydt[1], dydt[2]);
		printf("100%\n");
		exit(0);
#endif
	
	//Checked that this appears to behave correctly 2009-01-03
	
	//free things we've used, mostly all of the matrices
	gsl_matrix_free(Heff);
	gsl_matrix_free(HR);	
	gsl_matrix_free(Hani);
	gsl_matrix_free(HD);
	gsl_matrix_free(H);
	
	gsl_matrix_free(Hafc);
	
	gsl_matrix_free(p_mat);
	gsl_matrix_free(y_mat);
	gsl_matrix_free(Ht);

	return(GSL_SUCCESS);
}
