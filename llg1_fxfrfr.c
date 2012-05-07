/* 	llg1_fxfrfr.c - part of Spinsim_GSL Project
	This file only contains the function llg1()
	It is used by spinsim.c as the differential equation to solve for the magnetization
	of a 3 layer structure, with 1 fixed layer followed by 2 free layers.
	Written and engineered S. Moyerman
	University of California San Diego, ECE Department, Fullerton Lab
*/

//#define FXFRFR_HEFF_PRINT_AND_EXIT
//#define FXFRFR_SWAP_JEXCH_AND_HK

//Crude "scalar" of polarization; these are only used to scale an STT term entirely
//Leave these as 1.00!
#define POLARIZATION_LAYER_0 1.0000000
#define POLARIZATION_LAYER_1 1.0000000

/* llg1_fxfrfr() is a function that is the "differential equation" for the ODE solver
It takes parameters as defined by GSL odeiv system (http://www.gnu.org/software/gsl/manual/html_node/Defining-the-ODE-System.html)
t = the time position to use
y[] = the values of y at which our equation is being computed
dydt[] = the ar3lyr_print_heff_and_exitray in which to store the output
params = a pointer to a memory space of 23 pointers. See code below for how they are unpacked.
The function will return the values of dy/dt for each of the three dimensions STORED INTO dydt[]
The actual return value is GSL_SUCCESS for a valid execution */
int llg1_fxfrfr(double t, const double y[], double dydt[], void * params) {
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
	
	//ffprintf(stderr,stderr, "We got here 1\n");
	
	#ifdef FXFRFR_SWAP_JEXCH_AND_HK
	double* temp_af_coupling_js = (double*) malloc(6 * sizeof(double));
	for(count = 0; count < 6; count++) temp_af_coupling_js[count] = af_coupling_js[count];
	temp_af_coupling_js[0] = Hz[k];
	
	double Hz_k = 0; //We'll use this in place of Hz[k] for the effective field applied externally
	af_coupling_js = temp_af_coupling_js; //Use temp_af_coupling_js instead of the one provided (so we don't kill it for everyone else)
	#endif
	

//	Now we'll define and calculate variable's we'll actually use
//	Layer 0 only cares about polarization effect of layer 1 on it
	double Beta_div_P_and_j_10 = POLARIZATION_LAYER_0 * (GAMMA*HBAR)/(2*(dz[0]*1e-9)*MU0*Ms[0]*Q_E);
	double Beta = Beta_div_P_and_j_10*polarisationpp(y[0],y[1],y[2],y[3],y[4],y[5],pol[0],pol[1])*TimeDepFun(CALC_CURRENT, t, j[i], j_timedep_funs, t_relax);//remove me!
	double Beta_bar1 = Alpha[0]*Beta/( 1 + pow(Alpha[0],2) );
	double Beta_bar2 = Beta/(1+pow(Alpha[0],2));

	
//	Layer 1 cares about polarization effects from layer 0 and layer 2
//	Effect of layer 0 on layer 1
	double Beta_div_P_and_j_01 = POLARIZATION_LAYER_0 * (GAMMA*HBAR)/(2*(dz[1]*1e-9)*MU0*Ms[1]*Q_E);
	double Betap = -1 * Beta_div_P_and_j_01*polarisationpp(y[3],y[4],y[5],y[0],y[1],y[2],pol[1],pol[0])*TimeDepFun(CALC_CURRENT, t, j[i], j_timedep_funs, t_relax);
	double Beta_bar1p = Alpha[1]*Betap/( 1 + pow(Alpha[1],2) );
	double Beta_bar2p = Betap/(1+pow(Alpha[1],2));

//	Effect of layer 2 on layer 1
	double Beta_div_P_and_j_21 = POLARIZATION_LAYER_1 * (GAMMA*HBAR)/(2*(dz[1]*1e-9)*MU0*Ms[1]*Q_E);
	double Betap2 = Beta_div_P_and_j_21*polarisationpp(y[3],y[4],y[5],init_mag[6],init_mag[7],init_mag[8],pol[1],pol[2])*TimeDepFun(CALC_CURRENT, t, j[i], j_timedep_funs, t_relax);
	double Beta_bar1p2 = Alpha[1]*Betap2/( 1 + pow(Alpha[1],2) );
	double Beta_bar2p2 = Betap2/(1+pow(Alpha[1],2));

	/* Leave third layer fixed for now */ //IY comment wasn't closed

//	Re-define y as y_mat, a matrix with the values as different columns
	gsl_matrix* y_mat = gsl_matrix_calloc(3,1);
	gsl_matrix_set(y_mat, 0, 0, y[0]);
	gsl_matrix_set(y_mat, 1, 0, y[1]);
	gsl_matrix_set(y_mat, 2, 0, y[2]);
//	Re-define p as p_mat, a matrix with the values as different columns
//	This is our second "free" layer
	gsl_matrix* p_mat = gsl_matrix_calloc(3,1);
	gsl_matrix_set(p_mat, 0, 0, y[3]);
	gsl_matrix_set(p_mat, 1, 0, y[4]);
	gsl_matrix_set(p_mat, 2, 0, y[5]);
	
	// New for stationary layer
	// We grab the settings out of init_mag because this layer's magnetization will never change (given
	// that it is a fixed layer!)
	gsl_matrix* s_mat = gsl_matrix_calloc(3,1);
	gsl_matrix_set(s_mat, 0, 0, init_mag[6]);
	gsl_matrix_set(s_mat, 1, 0, init_mag[7]);
	gsl_matrix_set(s_mat, 2, 0, init_mag[8]);
	
	//ffprintf(stderr,stderr, "We got here 2\n");

	// With three layers, each one should have 3 demag contributions - one from each layer

//  Free Layer
//	Define H, the applied field, as a gsl_matrix, remember that it's (rows, columns)
	gsl_matrix* H = gsl_matrix_calloc(3,1);
	gsl_matrix_set(H, 0, 0, 0.0);
	gsl_matrix_set(H, 1, 0, 0.0);
	gsl_matrix_set(H, 2, 0, 0.0);
	
	if(sim_field_in_plane == 1) {
		#ifdef FXFRFR_SWAP_JEXCH_AND_HK
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz_k, hz_timedep_funs, t_relax));
		#else
		gsl_matrix_set(H, 0, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax));
		#endif
	}
	else {
		#ifdef FXFRFR_SWAP_JEXCH_AND_HK
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz_k, hz_timedep_funs, t_relax));
		#else
		gsl_matrix_set(H, 2, 0, TimeDepFun(CALC_FIELD, t, Hz[k], hz_timedep_funs, t_relax));
		#endif
	}
//	Anisotropic field
	gsl_matrix* Hani = gsl_matrix_calloc(3,1);
	gsl_matrix_set(Hani, 0, 0, H_K1[0]*y[0]*sin(thetaoffset));
	gsl_matrix_set(Hani, 1, 0, 0.0);
	gsl_matrix_set(Hani, 2, 0, H_K1[0]*y[2]*cos(thetaoffset)+H_K2[0]*y[2]*( pow(y[0],2) + pow(y[1],2) ) );
//	Self-demag field "Hani" for layer 0
	gsl_matrix* HD = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND[0], y_mat, 0, HD);
	gsl_matrix_scale(HD, (-1.0 * Ms[0])); //using anisotropy for layer 0
//	Mutual demag field "HR"
	gsl_matrix* HR = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR[0], p_mat, 0, HR);
	gsl_matrix_scale(HR, (-1.0 * Ms[1])); //using anisotropy for layer 1
//	(Anti-)Ferromagnetic Coupling Effective Field
//	Todo - there should be 1 such matrix for each of the 3 layer interactions!
	gsl_matrix* Hafc = gsl_matrix_calloc(3,1);
	//gsl_matrix_set(Hafc, 2, 0, af_coupling_js[0] * cos_angle_bw_two_vectors(y[0], y[1], y[2], init_mag[3], init_mag[4], init_mag[5])); //no good
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(Hafc, 0, 0, y[3] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 1, 0, y[4] * af_coupling_js[0] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 2, 0, y[5] * af_coupling_js[0] / (1e6 * Ms[0] * dz[0]*1e-9 * MU0));
	#else
	gsl_matrix_set(Hafc, 0, 0, y[3] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 1, 0, y[4] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc, 2, 0, y[5] * af_coupling_js[0] / (Ms[0] * dz[0]*1e-9));
	#endif
//	AF coupling from layer 2 ONTO layer 0
	gsl_matrix* Hafc2 = gsl_matrix_calloc(3,1);
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(Hafc2, 0, 0, init_mag[6] * af_coupling_js[1] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc2, 1, 0, init_mag[7] * af_coupling_js[1] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc2, 2, 0, init_mag[8] * af_coupling_js[1] / (1e6 * MU0 * Ms[0] * dz[0]*1e-9));
	#else
	gsl_matrix_set(Hafc2, 0, 0, init_mag[6] * af_coupling_js[1] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc2, 1, 0, init_mag[7] * af_coupling_js[1] / (Ms[0] * dz[0]*1e-9));
	gsl_matrix_set(Hafc2, 2, 0, init_mag[8] * af_coupling_js[1] / (Ms[0] * dz[0]*1e-9));
	#endif
	
//  Mutual demag field "HS"
//	Should this be between layer 0 and 2 or something else? I'm a little confused and should do
//	some thinking about this one.																			  
	gsl_matrix* HS = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR[1], s_mat, 0, HS);												   
	gsl_matrix_scale(HS, (-1.0 * Ms[2])); //using anisotropy for layer 2
	
	//How to create a random thermal "effective field" vector
	//Generate the random components
	double ht_comps0[3], ht_comps_len0;
	gsl_matrix* Ht0 = gsl_matrix_calloc(3,1);
	//int count;
	
	//Don't bother calculating if temperature is 0!
	if(temperature > 0.0) {
		for(count = 0; count < 3; count++) ht_comps0[count] = (1.0 * rand())/(1.0 * RAND_MAX);
		//Normalize said random components
		ht_comps_len0 = sqrt(pow(ht_comps0[0], 2) + pow(ht_comps0[1], 2) + pow(ht_comps0[2], 2));
		for(count = 0; count < 3; count++) ht_comps0[count] = ht_comps0[count] / ht_comps_len0;
		//Calculate the magnitude based on Boltzmann's constant etc, and multiply into components
		//For the dimension of the sample, we use dx * dy * (layer position of top layer + thickness of top layer)
		ht_comps_len0 = (temperature * K_BOLTZMANN) / (Ms[0] * (dx * 1e-9 * dy * 1e-9 * (dz[2]+layer_position_z[2]) * 1e-9));
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
	gsl_matrix_add(Heff, HS);
	gsl_matrix_add(Heff, Hafc); //commented out cause right now we don't want layer 0 to switch!
	gsl_matrix_add(Heff, Hafc2);
	gsl_matrix_add(Heff, Ht0);
	
	//ffprintf(stderr,stderr, "We got here 3\n");

//	  2nd free Layer (layer 1)
	gsl_matrix* HPani = gsl_matrix_calloc(3,1);
	gsl_matrix_set(HPani, 0, 0, H_K1[1]*y[5]*sin(thetaoffset));
	gsl_matrix_set(HPani, 1, 0, 0.0);
	gsl_matrix_set(HPani, 2, 0, H_K1[1]*y[5]*cos(thetaoffset) + H_K2[1]*y[5]*( pow(y[3],2) + pow(y[4],2) ));
	
	gsl_matrix* HPD = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND[1], p_mat, 0, HPD);
	gsl_matrix_scale(HPD, (-1.0 * Ms[1]));
	
	//Mutual demag tensor with layer 0 (layer 1 to layer 0)
	gsl_matrix* HPR = gsl_matrix_calloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR_down[0], y_mat, 0, HPR);
	gsl_matrix_scale(HPR, (-1.0 * Ms[0]));
		
	gsl_matrix* HPeff = gsl_matrix_calloc(3,1); //calloc sets initial values to 0
	gsl_matrix* HPS = gsl_matrix_calloc(3,1);
	
	//Mutual demag tensor with layer 2 (layer 1 to layer 2)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR[2], s_mat, 0, HPS);
	gsl_matrix_scale(HPS, (-1.0*Ms[2]));
	
//	(Anti-)Ferromagnetic Coupling Effective Field
//	Todo - there should be 1 such matrix for each of the 3 layer interactions!
	gsl_matrix* HPafc = gsl_matrix_calloc(3,1);
	
	//Fullerton comments that the AF coupling effective field should probably be J/M, where M = Ms * thickness_layer
	//There's an area term but since the area of both layers is equivalent, and we are doing a macro simulation,
	//we can probably ignore it.
	//We look at the Ms of the layer being acted upon, and the thickness of that layer
	
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(HPafc, 0, 0, y[0] * af_coupling_js[0] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 1, 0, y[1] * af_coupling_js[0] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 2, 0, y[2] * af_coupling_js[0] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	#else
	gsl_matrix_set(HPafc, 0, 0, y[0] * af_coupling_js[0] / (Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 1, 0, y[1] * af_coupling_js[0] / (Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc, 2, 0, y[2] * af_coupling_js[0] / (Ms[1] * dz[1]*1e-9));
	#endif
	
	gsl_matrix* HPafc2 = gsl_matrix_calloc(3,1);
	
	
	#ifdef JEXCHANGE_PEDANTIC
	gsl_matrix_set(HPafc2, 0, 0, init_mag[6] * af_coupling_js[3] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc2, 1, 0, init_mag[7] * af_coupling_js[3] / (1e6 * MU0 * Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc2, 2, 0, init_mag[8] * af_coupling_js[3] / (1e6 * Ms[1] * dz[1]*1e-9 * MU0));
	#else
	gsl_matrix_set(HPafc2, 0, 0, init_mag[6] * af_coupling_js[3] / (Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc2, 1, 0, init_mag[7] * af_coupling_js[3] / (Ms[1] * dz[1]*1e-9));
	gsl_matrix_set(HPafc2, 2, 0, init_mag[8] * af_coupling_js[3] / (Ms[1] * dz[1]*1e-9));
	#endif
	
	//How to create a random thermal "effective field" vector
	//Generate the random components
	double ht_comps1[3], ht_comps_len1;
	gsl_matrix* Ht1 = gsl_matrix_calloc(3,1);
	//int count;
	
	//Don't bother calculating if temperature is 0!
	if(temperature > 0.0) {
		for(count = 0; count < 3; count++) ht_comps1[count] = (1.0 * rand())/(1.0 * RAND_MAX);
		//Normalize said random components
		ht_comps_len1 = sqrt(pow(ht_comps1[0], 2) + pow(ht_comps1[1], 2) + pow(ht_comps1[2], 2));
		for(count = 0; count < 3; count++) ht_comps1[count] = ht_comps1[count] / ht_comps_len1;
		//Calculate the magnitude based on Boltzmann's constant etc, and multiply into components
		//For the dimension of the sample, we use dx * dy * (layer position of top layer + thickness of top layer)
		ht_comps_len1 = (temperature * K_BOLTZMANN) / (Ms[1] * (dx * 1e-9 * dy * 1e-9 * (dz[2]+layer_position_z[2]) * 1e-9));
		for(count = 0; count < 3; count++) ht_comps1[count] = ht_comps1[count] * ht_comps_len1;
		//Put all of random components into a matrix
		gsl_matrix_set(Ht1, 0, 0, ht_comps1[0]);
		gsl_matrix_set(Ht1, 1, 0, ht_comps1[1]);
		gsl_matrix_set(Ht1, 2, 0, ht_comps1[2]);
	}
	
	gsl_matrix_add(HPeff, H);
	gsl_matrix_add(HPeff, HPani);
	gsl_matrix_add(HPeff, HPR);
	gsl_matrix_add(HPeff, HPD);
	gsl_matrix_add(HPeff, HPS);
	gsl_matrix_add(HPeff, HPafc);
	gsl_matrix_add(HPeff, HPafc2);
	gsl_matrix_add(HPeff, Ht1);

#ifdef FXFRFR_HEFF_PRINT_AND_EXIT
	if(t >= 2e-8) {
	fprintf(stderr,"Heff = [%.2e, %.2e, %.2e]\n", gsl_matrix_get(Heff,0,0), gsl_matrix_get(Heff,1,0), gsl_matrix_get(Heff,2,0));
//	fprintf(stderr,"Beta_div_P_and_j = %.2e, polarisation = %.2e\n", Beta_div_P_and_j, polarisation(y[0],y[1],y[2],y[3],y[4],y[5]));

	fprintf(stderr,"Hani = [%.2e, %.2e, %.2e]\nHD = [%.2e, %.2e, %.2e]\nHR = [%.2e, %.2e, %.2e]\nHS = [%.2e, %.2e, %.2e]\n", 
		gsl_matrix_get(Hani, 0, 0), gsl_matrix_get(Hani, 1, 0), gsl_matrix_get(Hani, 2, 0), 
		gsl_matrix_get(HD, 0, 0), gsl_matrix_get(HD, 1, 0), gsl_matrix_get(HD, 2, 0), 
		gsl_matrix_get(HR, 0, 0), gsl_matrix_get(HR, 1, 0), gsl_matrix_get(HR, 2, 0), 
		gsl_matrix_get(HS, 0, 0), gsl_matrix_get(HS, 1, 0), gsl_matrix_get(HS, 2, 0));
		
	/*fprintf(stderr,"HR = [%.2e, %.2e, %.2e]\nHS = [%.2e, %.2e, %.2e]\nHPR = [%.2e, %.2e, %.2e]\nHPS = [%.2e, %.2e, %.2e]\n", 
		gsl_matrix_get(HR, 0, 0), gsl_matrix_get(HR, 1, 0), gsl_matrix_get(HR, 2, 0), 
		gsl_matrix_get(HS, 0, 0), gsl_matrix_get(HS, 1, 0), gsl_matrix_get(HS, 2, 0), 
		gsl_matrix_get(HPR, 0, 0), gsl_matrix_get(HPR, 1, 0), gsl_matrix_get(HPR, 2, 0), 
		gsl_matrix_get(HPS, 0, 0), gsl_matrix_get(HPS, 1, 0), gsl_matrix_get(HPS, 2, 0));*/
		
	fprintf(stderr,"Hafc = [%.2e, %.2e, %.2e]\nHafc2 = [%.2e, %.2e, %.2e]\n", 
		gsl_matrix_get(Hafc, 0, 0), gsl_matrix_get(Hafc, 1, 0), gsl_matrix_get(Hafc, 2, 0), 
		gsl_matrix_get(Hafc2, 0, 0), gsl_matrix_get(Hafc2, 1, 0), gsl_matrix_get(Hafc2, 2, 0));
	
	fprintf(stderr,"HPeff = [%.2e, %.2e, %.2e]\n", gsl_matrix_get(HPeff,0,0), gsl_matrix_get(HPeff,1,0), gsl_matrix_get(HPeff,2,0));
//	fprintf(stderr,"Beta_div_P_and_j = %.2e, polarisation = %.2e\n", Beta_div_P_and_j, polarisation(y[0],y[1],y[2],y[3],y[4],y[5]));

	fprintf(stderr,"HPani = [%.2e, %.2e, %.2e]\nHPD = [%.2e, %.2e, %.2e]\nHPR = [%.2e, %.2e, %.2e]\nHPS = [%.2e, %.2e, %.2e]\n", 
		gsl_matrix_get(HPani, 0, 0), gsl_matrix_get(HPani, 1, 0), gsl_matrix_get(HPani, 2, 0), 
		gsl_matrix_get(HPD, 0, 0), gsl_matrix_get(HPD, 1, 0), gsl_matrix_get(HPD, 2, 0), 
		gsl_matrix_get(HPR, 0, 0), gsl_matrix_get(HPR, 1, 0), gsl_matrix_get(HPR, 2, 0), 
		gsl_matrix_get(HPS, 0, 0), gsl_matrix_get(HPS, 1, 0), gsl_matrix_get(HPS, 2, 0));
		
	/*fprintf(stderr,"HR = [%.2e, %.2e, %.2e]\nHS = [%.2e, %.2e, %.2e]\nHPR = [%.2e, %.2e, %.2e]\nHPS = [%.2e, %.2e, %.2e]\n", 
		gsl_matrix_get(HR, 0, 0), gsl_matrix_get(HR, 1, 0), gsl_matrix_get(HR, 2, 0), 
		gsl_matrix_get(HS, 0, 0), gsl_matrix_get(HS, 1, 0), gsl_matrix_get(HS, 2, 0), 
		gsl_matrix_get(HPR, 0, 0), gsl_matrix_get(HPR, 1, 0), gsl_matrix_get(HPR, 2, 0), 
		gsl_matrix_get(HPS, 0, 0), gsl_matrix_get(HPS, 1, 0), gsl_matrix_get(HPS, 2, 0));*/
		
	fprintf(stderr,"HPafc = [%.2e, %.2e, %.2e]\nHPafc2 = [%.2e, %.2e, %.2e]\n", 
		gsl_matrix_get(HPafc, 0, 0), gsl_matrix_get(HPafc, 1, 0), gsl_matrix_get(HPafc, 2, 0), 
		gsl_matrix_get(HPafc2, 0, 0), gsl_matrix_get(HPafc2, 1, 0), gsl_matrix_get(HPafc2, 2, 0));
		
	fprintf(stderr,"t = %.2e\nL1 = (%.2e, %.2e, %.2e)\nL2 = (%.2e, %.2e, %.2e)\n",t,y[0],y[1],y[2],y[3],y[4],y[5]);
	printf("100\n");
	exit(0);
	}
#endif

//	In fact, one has to integrate three differential equations,
//	one for each component of m

	// Free Layer Dynamics
	dydt[0] = -1.0 * Gamma_bar[0] * (y[1]*gsl_matrix_get(Heff, 2, 0)-y[2]*gsl_matrix_get(Heff, 1, 0)) - Alpha_bar[0]*((y[0]*y[0]-1)*gsl_matrix_get(Heff, 0, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 1, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 2, 0)) + Beta_bar1 * (y[1]*y[5] - y[2]*y[4]) - Beta_bar2 * ((y[0]*y[0]-1)*y[3]+y[0]*y[1]*y[4]+y[0]*y[2]*y[5]);
	dydt[1] = -1.0 * Gamma_bar[0]*(y[2]*gsl_matrix_get(Heff, 0, 0)-y[0]*gsl_matrix_get(Heff, 2, 0))-Alpha_bar[0]*((y[1]*y[1]-1)*gsl_matrix_get(Heff, 1, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 2, 0))+Beta_bar1*(y[2]*y[3]-y[0]*y[5]) - Beta_bar2*((y[1]*y[1]-1)*y[4]+y[0]*y[1]*y[3]+y[1]*y[2]*y[5]);
	dydt[2] = -1.0 * Gamma_bar[0]*(y[0]*gsl_matrix_get(Heff, 1, 0)-y[1]*gsl_matrix_get(Heff, 0, 0))-Alpha_bar[0]*((y[2]*y[2]-1)*gsl_matrix_get(Heff, 2, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 1, 0))+Beta_bar1*(y[0]*y[4]-y[1]*y[3])-Beta_bar2*((y[2]*y[2]-1)*y[5]+y[0]*y[2]*y[3]+y[1]*y[2]*y[4]);
	
	dydt[3] = -1.0 * Gamma_bar[1]*(y[4]*gsl_matrix_get(HPeff, 2, 0)-y[5]*gsl_matrix_get(HPeff, 1, 0)) - Alpha_bar[1]*((y[3]*y[3]-1)*gsl_matrix_get(HPeff, 0, 0)+y[3]*y[4]*gsl_matrix_get(HPeff, 1, 0)+y[3]*y[5]*gsl_matrix_get(HPeff, 2, 0))+Beta_bar1p*(y[4]*y[2]-y[5]*y[1])-Beta_bar2p*((y[3]*y[3]-1)*y[0]+y[3]*y[4]*y[1]+y[3]*y[5]*y[2]) + Beta_bar1p2*(y[4]*init_mag[8]-y[5]*init_mag[7])-Beta_bar2p2*((y[3]*y[3]-1)*init_mag[6]+y[3]*y[4]*init_mag[7]+y[3]*y[5]*init_mag[8]);
	dydt[4] = -1.0*Gamma_bar[1]*(y[5]*gsl_matrix_get(HPeff, 0, 0)-y[3]*gsl_matrix_get(HPeff, 2, 0))-Alpha_bar[1]*((y[4]*y[4]-1)*gsl_matrix_get(HPeff, 1, 0)+y[3]*y[4]*gsl_matrix_get(HPeff, 0, 0)+y[4]*y[5]*gsl_matrix_get(HPeff, 2, 0))+Beta_bar1p*(y[5]*y[0]-y[3]*y[2])-Beta_bar2p*((y[4]*y[4]-1)*y[1]+y[3]*y[4]*y[0]+y[4]*y[5]*y[2])+Beta_bar1p2*(y[5]*init_mag[6]-y[3]*init_mag[8])-Beta_bar2p2*((y[4]*y[4]-1)*init_mag[7]+y[3]*y[4]*init_mag[6]+y[4]*y[5]*init_mag[8]);
	
	dydt[5] = -1.0*Gamma_bar[1]*(y[3]*gsl_matrix_get(HPeff, 1, 0)-y[4]*gsl_matrix_get(HPeff, 0, 0))-Alpha_bar[1]*((y[5]*y[5]-1)*gsl_matrix_get(HPeff, 2, 0)+y[3]*y[5]*gsl_matrix_get(HPeff, 0, 0)+y[4]*y[5]*gsl_matrix_get(HPeff, 1, 0))+Beta_bar1p*(y[3]*y[1]-y[4]*y[0])-Beta_bar2p*((y[5]*y[5]-1)*y[2]+y[3]*y[5]*y[0]+y[4]*y[5]*y[1])+Beta_bar1p2*(y[3]*init_mag[7]-y[4]*init_mag[6])-Beta_bar2p2*((y[5]*y[5]-1)*init_mag[8]+y[3]*y[5]*init_mag[6]+y[4]*y[5]*init_mag[7]);
	
	//fprintf(stderr,"dydy[] = %.3e %.3e %.3e %.3e %.3e %.3e\n", dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5]);
	//exit(0);
	//fprintf(stderr,"Gamma bar = %f; y[] = %.3e %.3e %.3e %.3e %.3e %.3e;\nHeff = %.3e %.3e %.3e; Alpha bar = %.3e; Beta_bar =- %.3e %.3e; Beta = %f / %f\npol = %f, j = %f, i = %d\n", 
		//Gamma_bar, y[0], y[1], y[2], y[3], y[4], y[5], gsl_matrix_get(Heff, 0, 0), gsl_matrix_get(Heff, 1, 0), gsl_matrix_get(Heff, 2, 0),
		//Alpha_bar, Beta_bar1, Beta_bar2, Beta, Beta_div_P_and_j, polarisation(y[0],y[1],y[2],y[3],y[4],y[5]), j[i], i);
	
	//exit(0);
	
	//Checked that this appears to behave correctly 2009-01-03
	
	//ffprintf(stderr,stderr, "We got here 4\n");
	
	//free things we've used, mostly all of the matrices
	gsl_matrix_free(Heff);
	gsl_matrix_free(HR);	
	gsl_matrix_free(Hani);
	gsl_matrix_free(HD);
	gsl_matrix_free(H);
	gsl_matrix_free(HS);
	gsl_matrix_free(p_mat);
	gsl_matrix_free(y_mat);
	gsl_matrix_free(HPeff);
	gsl_matrix_free(HPD);
	gsl_matrix_free(HPani);
	gsl_matrix_free(HPR);
	gsl_matrix_free(HPS);
	gsl_matrix_free(s_mat);
	
	gsl_matrix_free(Hafc);
	gsl_matrix_free(HPafc2);
	gsl_matrix_free(HPafc);
	gsl_matrix_free(Hafc2);
	
	gsl_matrix_free(Ht0);
	gsl_matrix_free(Ht1);
	
	#ifdef FXFRFR_SWAP_JEXCH_AND_HK
	free(temp_af_coupling_js);
	#endif
	
	//(stderr, "We got here 5\n");

	return(GSL_SUCCESS);
}
