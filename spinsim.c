
/* 	spinsim.c
	This file is the top level code for the spinsim simulator. It allows simulation of the behavior of a free layer
	across a grid of externally applied current/H-field conditions. The function to call is Hsweep()
	For greater detail, see documentation PDF for this project.
	Note: This was written and tested on Linux 2.6.18 x86-64 with GSL version 1.12, it is recommended to use
	the same versions to prevent compatibility issues.
	Note: If compiling using gcc, be sure to use the '-pthread -lgsl -lgslcblas -lm' options on the command line.
	Written I. Yulaev 2010-01-xx (iyulaev@ucsd.edu)
	Engineered S. Moyerman (smoyerman@ucsd.edu)
*/

//SPECIAL, USER-SPECIFIED DEFINES BELOW:

//Define the following if you want to calculate GMR, rather than just free layer magnetization
//#define CALCULATE_GMR 0

//Define the following if you want to average the magnetization over time, rather than just take the final magnetization
//#define AVERAGE_MAGNETIZATION 1

//Set this to 1 if you want to look at layer 2's magnetization instead of layer 1's
//This also changes what layer gets its magnetization inverted when we run sweep-up vs. sweep-down
#define LOOK_AT_LAYER_2 0

//Un-comment (enable) this define if you want both layers 0 and 1 to have their initial magnetizations
//switched (flipped across x-y plane) when running 'sweep-up' simulation
#define SWITCH_BOTH_LAYERS

//Define this to have spinsim spit out when sweep-up and sweep-down magnetization (or GMR, if CALCULATE_GMR defined)
//switches from being positive to negative or visa versa
//Also it'll spit out the middle of the switching window
//Warning, experimental, may not work well
//Only use when number of threads is 1
//#define TRACK_SWITCHING_POINT 1

//The below defines allow for the tracking of tsettling. This is useful when doing measurements
//of alpha within the system since alpha ~ 1/tsettling
//#define TRACK_TSETTLING_X
//#define TSETTLING_CRITERIA 0.1
//#define TRACK_TSETTLING_Z
//#define TSETTLING_CRITERIA_Z -0.95

//#define POLARISATION_XIAO
#define POLARISATION_SLONCZEWSKI

//Uncomment this define if you want to re-start the y-values from the 'initial values' 
//instead of sweeping across current each time
#define RESTART_Y_FROM_INIT_VALS
//#define SWEEP_CURRENT
//#define SWEEP_FIELD_AND_CURRENT

//Uncomment the below to use uJ/m**2 as the units for exchange coupling.
//If commented, "legacy" units will be used
#define JEXCHANGE_PEDANTIC

//END SPECIAL USER-SPECIFIED DEFINES


//Defines for some physical constants
#define K_BOLTZMANN 1.38065e-23
#define GAMMA 2.21e5 //Gyromagnetic ratio, 2.2127616 x 1e5 s^(-1).A^(-1).m
#define MU0 1.3e-6 //Magnetic permeability of vaccuum, H/m, 4*PI*1e-7 H/m
#define Q_E 1.602e-19 //Electron charge, -1.602176 x 1e-19 C
#define HBAR 1.054e-34 //Reduced Planck constant, 1.054571 x 1e-34 J.s

//Defines for integrator parameters
#define INTEGRATOR_TO_USE gsl_odeiv_step_rk4
//Defines for integrator precision
#define ABS_PREC 1e-6			//Absolute precision for integrator to aim for
#define REL_PREC 1e-6		//Relative precision for integrator to aim for
#define PREC_MODIFIER 100.0		//For some simulation types (i.e. Fixed - Free - Free or upside-down of that) the integrator will never
								//finish with certain precision values. We modify the precision value by a scalar in those cases.
#define INIT_STEP_SZ 100000.0	//Define the proportion that the initial step will be of t_f, i.e. it will be t_f/INIT_STEP_SZ
								//This doesn't really make a big difference...
#define INT_RENORM_VECTOR		//Define whether or not to re-normalize the vector after every integration
#define INT_RENORM_VECTOR_TOL 0.9999999 //Define the tolerance for re-normalizing the vector
#define INT_MIN_STEP_SZ	1e-18	//Define the minimum step size, in seconds

//Internal defines, DON'T TOUCH! Unless you're modifying program infrastructure
#define NUMARGS 60
#define NOISINESS 1		//How talkative? 1 is "as loud as possible", a large number makes it quiet (like 10000)

//Debugging output
#define INT_TIME_DEBUG 1 //print out debug messages if the integrator takes a while for a particular step
#define INT_TIME_SECONDS 30.0 //takes at least this long before debug message prints out
//#define MIN_VEC_LEN_PRINT	0.999 //below this value, the integrator will print the current vector length to STDERR

#include <stdio.h>
#include <math.h>
#include <time.h>

//Includes from GSL Library
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_odeiv.h>

//Includes for PThreads API
#include <sys/types.h>
#include <pthread.h>

//Include helper functions like demagtensor and polarisation
#include "spin_helpers.c"
//Include the functions for calculating time dependent field and current
#include "timedepfuns.c"
//Include the differential equation functions
#include "llg1.c"
#include "llg1_6.c"
#include "llg1_fxfrfr.c"

/* Hsweep() calculates the magnetization of a free layer across a 2D plane of external field & applied current conditions
Parameters:
(1) Theta offset (offset, in degrees, that the applied H-field is from the Z axis of the sample
(2) Percent that K2 (2nd order anisotropy) is of K1
(3) Alpha (gilbert's damping constant)
(4) (double*) of size 4, Layer1 through layer 4 saturation magnetization
(5) (double*) of size 4, Layer 1 through layer 4 anisotropy (K1)
(6) (double*) of size 4, layer position on z axis for layers 1-4
(7) Simulation type, lower 8 bits used. Each 2-bit pair encodes layers 1-4 (4 - LS bits) types. 00 means 
	disabled, 01 means fixed, 10 means free
(8-10) Minimum current density, maximum current density, number of current density points to plot against (units are A/m**2)
(11-13) Minimum applied field, maximum applied field, number of applied field points to plot against (units are teslas)
(14) Integration time, in seconds
(15) J time dependent functions, bit-wise OR'd (see documentation)
(16) H time dependent functions, bit-wise OR'd (see documentation)
(17-18) Sample size, x and y axis
(19) Temperature of the sample, for "thermal noise" - use 0 to ignore temperature-dependent noise
(20) (double*) of size 12, the initial magnetizations (x,y,z) triples for layers 1-4
(21) (double*) of size 4, the thicknesses for layers 1-4
(22) Pointer for j (output)
(23) Pointer for hz (output)
(24) Pointer for magsweepup data (output)
(25) Pointer for magsweepdown data (output)
(26) Number of threads to use
(27) The number of this thread
(28) The AF coupling constants (doubles) for layers 1-2, 1-3, 1-4, 2-3, 2-4, 3-4

Remarks: If n == 1 and n_field==1 then it saves every time step for the single point solution 
and returns a pointer to the memory space where the result is stored (dynamically allocated)

Return value: Returns NULL, or if called with n = 1, n_field = 1, it returns a pointer to the 
memory space where the output data is stored. 
*/
double* Hsweep(double thetaoff, double percK2, double * Alpha, double* Ms, double* K1, double* layer_position_z,
unsigned int symtype, double j_min, double j_max, int n, double hz_min, double hz_max, int n_field, double t_f,
unsigned int j_timedep_funs, unsigned int hz_timedep_funs, double dx, double dy, double temperature, 
double* init_mag, double* dz, double* j, double* Hz, double* magsweepup, double* magsweepdown, int numthreads,
int thread_no, double* af_coupling_js, double* pol, double t_relax){
	
//	save start and stop times to benchmark the program speed
	time_t start, stop, printed_last, temp;
	time(&start); 
	
//	number of dimensions, we'll set this later
	int numdims;
	
//	temporary loop variables
	int count, count2;
	
//	Save whether or not to simulate in plane (and remove it from simtype)
	int sim_field_in_plane = (symtype >> 31)&1;
	if(sim_field_in_plane == 1) {
		symtype -= (1 << 31);
		fprintf(stderr, "GOING TO RUN SIM FIELD IN PLANE!\n");
	}
	
//	These are used to calculate the settingling time function
	
	#ifdef TRACK_TSETTLING_X
	double last_three_lyr1_max[3];
	double last_three_lyr1_pts[3];
	double last_three_lyr2_max[3];
	double last_three_lyr2_pts[3];
	double tsettling_lyr1;
	double tsettling_lyr2;
	#endif
	#ifdef TRACK_TSETTLING_Z
	double tsettling_lyr1;
	double tsettling_lyr2;
	#endif
	
//	make all initial magnetization vectors unit-sized
	for(count = 0; count < 4; count++) {
		if( ((symtype >> (2*count)) & 3) != 0 ) {
			double vec_size = sqrt( pow(( init_mag[3*count] ),2) + pow(( init_mag[3*count+1] ),2) + pow(( init_mag[3*count+2] ),2) );
			init_mag[3*count] = init_mag[3*count] / vec_size;
			init_mag[3*count+1] = init_mag[3*count+1] / vec_size;
			init_mag[3*count+2] = init_mag[3*count+2] / vec_size;
		}
	}

//	Integration times
	double t_0 = 0;
//	Precision for integrator
	double absolute_precision = ABS_PREC * ((symtype == 26 || symtype == 41)?PREC_MODIFIER : 1.0);
	double relative_precision = REL_PREC * ((symtype == 26 || symtype == 41)?PREC_MODIFIER : 1.0);
	
//	Return value (double*) in case we're doing a single point run, with n == 1 and n_field == 1
	double * retval = malloc(100000*4*sizeof(double)); //we set it to 100,000 doubles arbitrarily
//	We'll re-alloc in the integration loop as necessary, so size is somewhat arbitrary
	int retval_size = 100000;
	int num_evolutions = 0;
	
//	These can potentially be used to track the switching point, if TRACK_SWITCHING_POINT is defined
	#ifdef TRACK_SWITCHING_POINT
	double su_switch_hz, su_switch_j, sd_switch_hz, sd_switch_j;
	#endif
	
//	Self-Demag Tensors, make one for every layer 1-4
	gsl_matrix * ND[4];
	for(count = 0; count < 4; count++) ND[count] = NULL;
	for(count = 0; count < 4; count++) {
//		check to make sure we're actually to use that layer!
		if( ((symtype >> (2*count)) & 3) != 0 ) {
			ND[count] = gsl_matrix_alloc(3,3);
			
//			demagtensor takes parameters X,Y,Z which is layer position, and dx dy dz, which is layer extents
			demagtensor(0.0, 0.0, 0.0, dx, dy, dz[count], ND[count]);
		}
	}

//	Mutual Demagnetization Tensors
//	One will be intialized for every PAIR of used layers
//	First, we create a tensor for every layer n to every layer n-1
	gsl_matrix * NR[6];
	for(count = 0; count < 6; count++) NR[count] = NULL;
	int mat_count = 0;
	for(count = 0; count < 3; count++) {
		for(count2 = count+1; count2 < 4; count2++) {
//			Check to make sure that the pair of layers is actually used!
			if( (( (symtype >> (2*count)) & 3) != 0 ) && ( ((symtype >> (2*count2)) & 3) != 0 ) ) {
				NR[mat_count] = gsl_matrix_calloc(3,3);

//				Calculate demagtensor between layers count2 (as the layer for the position) and count (as the layer for the extents)
//				TODO: Ivan figure out if this is reasonable
				demagtensor(0.0,0.0,abs(layer_position_z[count2]-layer_position_z[count]), dx, dy, dz[count], NR[mat_count++]);
			}
		}
	}
	
	//	First, we create a tensor for every layer n to every layer n+1
	gsl_matrix * NR_down[6];
	for(count = 0; count < 6; count++) NR_down[count] = NULL;
	mat_count = 0;
	for(count = 0; count < 3; count++) {
		for(count2 = count+1; count2 < 4; count2++) {
//			Check to make sure that the pair of layers is actually used!
			if( (( (symtype >> (2*count)) & 3) != 0 ) && (( (symtype >> (2*count2)) & 3) != 0 ) ) {
				NR_down[mat_count] = gsl_matrix_calloc(3,3);

//				Calculate demagtensor between layers count (as the layer for the position) and count2 (as the layer for the extents)
//				TODO: Ivan figure out if this is reasonable
//				Shoudl it be the difference between the two layer positions, maybe?
				demagtensor(0.0,0.0,abs(layer_position_z[count]-layer_position_z[count2]), dx, dy, dz[count2], NR_down[mat_count++]);
			}
		}
	}
	
//	gsl_matrix * NR = gsl_matrix_alloc(3,3);
//	demagtensor(X,Y,Z,dx,dy,dz,NR);
//	Note: For infinite lateral size use NR_xx = 0; NR_yy = 0; NR_zz = 0;

//	Self Demagnetization Matrix
//	X = 0; Y = 0; Z = 0;
//	gsl_matrix * ND = gsl_matrix_alloc(3,3);
//	demagtensor(X,Y,Z,dx,dy,dz,ND);
//	Note: For infinite lateral size use ND_xx = 0; ND_yy = 0; ND_zz = 1;

//	Define the normal LLG parameters for the problem
	double Gamma = 2.21e5; 			//Gyromagnetic ratio, 2.2127616 x 1e5 s^(-1).A^(-1).m
	double mu0 = 4*PI*1e-7; 		//Magnetic permeability of vaccuum, H/m
	double e = 1.602e-19; 			//Electron charge, -1.602176 x 1e-19 C
	double hbar = 1.054e-34; 		//Reduced Planck constant, 1.054571 x 1e-34 J.s
	
//	Generate K2 based on K1[] and count
	double K2[4]; //Free layer anisotropy constant
	for(count = 0; count < 4; count++) K2[count] = (percK2/100)*K1[count];

//	Other constants
		
//	Slonczewski's spin transfert torque coefficient Beta=P*j*Beta_div_P_and_j
//	Ivan: This is a hack for 2 layer problem only
//	Don't forget to convert dz into meters, as it's given in nanometers (I hate this!)
	double Beta_div_P_and_j = (Gamma*hbar)/(2*(dz[0]*1e-9)*mu0*Ms[0]*e);
	
	//H_K2 = 0; //This doesn't feel right, but it was originally uncommented
	double Gamma_bar[4];
	for(count=0; count < 4; count++) Gamma_bar[count] = Gamma/(1+pow(Alpha[count],2) );
	double Alpha_bar[4];
	for(count=0; count < 4; count++) Alpha_bar[count] = Alpha[count]*Gamma/( 1 + pow(Alpha[count],2) );

//	Some useful parameters
	double H_K1[4]; //coefficient for field due to 1st order anisotropy
	double H_K2[4]; //coefficient for field due to 2nd order anisotropy
	double norm[4]; //normalized 1st order field due to anisotropy term

	for(count = 0; count < 4; count++) {
		H_K1[count] = 2*K1[count]/mu0/Ms[count];				//Free layer anisotropy field, A/m
		H_K2[count] = 4*K2[count]/mu0/Ms[count];				//Free layer anisotropy field, A/m
		
		norm[count] = H_K1[count] / sqrt( pow(H_K1[count],2)+pow(H_K2[count],2) );
		H_K1[count] = norm[count]*H_K1[count];
		H_K2[count] = norm[count]*H_K2[count];
	}
	
//	stuff for the 3rd layer
//	Define params for the reference layer - for reasonable dynamics, make 1.5x K1
//	Ivan: No such thing anymore! We don't differntial between pinned and unpinned layers! It's all done in
//	the parameters now

//	double Kp1 = 3e5*10;          // pinned layer anisotropy
//	double HP_K1 = 2*Kp1/mu0/Ms;                        //Pinned layer anisotropy field, A/m

//	Initialize our applied field vector
//	Applied field will go from (-HZ_MAX) to (HZ_MAX) with n_field steps in between
	make_linspace(hz_min,hz_max,n_field,Hz);
	
//	Current Density Vector
//	Applied current density will go from (-J_MAX) to (J_MAX) with n_field steps in between
	make_linspace(j_min,j_max,n,j); 

//	Offset for easy axis of free FM layer, convert to radians
	double thetaoffset = thetaoff*PI/180;
	
//	Loop variables
	int i, k;
	char percent[2] = "%";
	double curr_value=0.0;
	double y_res_5, y_res_2;
	#ifdef TRACK_SWITCHING_POINT
	double* prev_value_sd = (double*) calloc(n_field * n * sizeof(double),1);
	double* prev_value_su =  (double*) calloc(n_field * n * sizeof(double),1);
	int* su_init = (int*) calloc (n_field* n *sizeof(int),1);
	int* sd_init = (int*) calloc (n_field* n *sizeof(int),1);
	#endif
	
//	We'll use this to gauge when the last time we gave a status update was
	time(&printed_last);
	
	#ifdef SWEEP_FIELD_AND_CURRENT
//	Make space for the output of the integration
	double y_result[6];
	double y_result_up[6];
	double y_result_down[6];
//	Set initial condition (we only set this once, at the beginning of the integration
	y_result_down[0] = init_mag[0];
	y_result_down[1] = init_mag[1];
	y_result_down[2] = init_mag[2];
	y_result_down[3] = init_mag[3];
	y_result_down[4] = init_mag[4];
	y_result_down[5] = init_mag[5];
	
	//		Set initial condition (multiply one of them by -1, depending on what the "free" layer we're looking at is
	#if LOOK_AT_LAYER_2 == 0
	y_result_up[0] = init_mag[0];
	y_result_up[1] = init_mag[1];
	y_result_up[2] = (-1.0) * init_mag[2];
	y_result_up[3] = init_mag[3];
	y_result_up[4] = init_mag[4];
		#ifdef SWITCH_BOTH_LAYERS
			y_result_up[5] = (-1.0) * init_mag[5];
		#endif
		#ifndef SWITCH_BOTH_LAYERS
			y_result_up[5] = init_mag[5];
		#endif
	#elif LOOK_AT_LAYER_2 == 1
	y_result_up[0] = init_mag[0];
	y_result_up[1] = init_mag[1];
		#ifdef SWITCH_BOTH_LAYERS
			y_result_up[2] = (-1.0) * init_mag[2];
		#endif
		#ifndef SWITCH_BOTH_LAYERS
			y_result_up[2] = init_mag[2];
		#endif
	y_result_up[3] = init_mag[3];
	y_result_up[4] = init_mag[4];
	y_result_up[5] = (-1.0) * init_mag[5];
	#endif
	 
	#endif
	
	for(i = thread_no; i < n; i+=numthreads) {
		time(&temp);
		if((thread_no == 0) && ( (((int)temp) - ((int)printed_last)) >= NOISINESS) ) {
			printf("%.1f%s\n", ((float) i) / ((float) n) * 100.0, percent);
			fflush(stdout);
			time(&printed_last);
		}
		
		#ifdef SWEEP_CURRENT
//		Make space for the output of the integration
		double y_result[6];
//		Set initial condition (we only set this once, at the beginning of the integration
		y_result[0] = init_mag[0];
		y_result[1] = init_mag[1];
		y_result[2] = init_mag[2];
		y_result[3] = init_mag[3];
		y_result[4] = init_mag[4];
		y_result[5] = init_mag[5]; 
		#endif
		
//		sweep down
		for(k = (n_field-1); k >= 0; k--) {
		//Ivan: Stopped here 2010-02-14, need to figure out which parameter can be condensed
//			put together a parameter list... (this whole thing needs to be re-written)
			double ** param_list = (double **) malloc(27 * sizeof(double *));
			param_list[0] = init_mag;
			param_list[1] = Hz;
			param_list[2] = j;
			param_list[3] = &Beta_div_P_and_j;
			param_list[4] = Alpha;
			param_list[5] = Alpha_bar;
			param_list[6] = H_K1;
			param_list[7] = H_K2;
			param_list[8] = &thetaoffset;
			param_list[9] = (double *) ND;
			param_list[10] = (double *) NR;
			param_list[11] = (double *) NR_down;
			param_list[12] = Ms;
			param_list[13] = Gamma_bar;
			param_list[14] = (double *) &k; //yes, k is an int, so what
			param_list[15] = (double *) &i; //yes, i is an int, this may not be portable but should be fine for x86(-64)
			param_list[16] = (double*) &j_timedep_funs;
			param_list[17] = (double*) &hz_timedep_funs;
			param_list[18] = (double*) &temperature;
			param_list[19] = &dx;
			param_list[20] = &dy;
			param_list[21] = dz;
			param_list[22] = layer_position_z;
			param_list[23] = af_coupling_js;
			param_list[24] = pol;
			param_list[25] = &t_relax;
			param_list[26] = (double*)&sim_field_in_plane;
			
//			set up the ode system (most of this code is lifted from the GSL documentation and just re-written)
//			http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html
//			we don't give it a jacobian because we're going to use runge-kutta algorithm which doesn't need it

			gsl_odeiv_system* sys_ptr;
			
			//2D simulation case
			if(symtype == 6 || symtype == 9) {
				gsl_odeiv_system sys = {llg1, NULL, 3, param_list};
				numdims = 3;
				sys_ptr = &sys;
			}
			else if(symtype == 10) {
				gsl_odeiv_system sys = {llg1_6, NULL, 6, param_list};
				numdims = 6;
				sys_ptr = &sys;
			}
			//3D, FxFrFr or FrFrFx
			else if(symtype == 41 || symtype == 26) {
				gsl_odeiv_system sys = {llg1_fxfrfr, NULL, 6, param_list};
				numdims = 6;
				sys_ptr = &sys;
			}
			//If we're doing a simulation with the fixed layer on the BOTTOM reverse the direction of interactions
			//by reversing the down->up and up->down mutual demag tensors
			if(symtype == 9 || symtype == 41) {
				double* temp = param_list[10];
				param_list[10] = param_list[11];
				param_list[11] = param_list[10];
			}
			
//			define step type
//			So far, gsl_odeiv_step_rk4 seems to be the best, but gsl_odeiv_step_rk8pd, gsl_odeiv_step_rk2 isn't bad either. R-K (4,5) appears to suck.
			const gsl_odeiv_step_type * T = INTEGRATOR_TO_USE;
//			define the "stepper", give it the step type and # of dimensions
			gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, numdims);
//			Define controller
//			Ivan: Stephanie, please examine this - do we need to use the _y or _yp version? Also parameters to this controller, I just guessed.
//			see http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html
			//gsl_odeiv_control * c = gsl_odeiv_control_y_new (absolute_precision, relative_precision);
			gsl_odeiv_control * c = gsl_odeiv_control_standard_new (absolute_precision, relative_precision, 1.0, 1.0);
//			Allocate the function we're going to use for evolution
			gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (numdims);
//			Define the starting time for our integration
			double t = t_0;
//			Define the initial step size
			double h = t_f / INIT_STEP_SZ;
			
			#ifdef RESTART_Y_FROM_INIT_VALS
//			Make space for the output of the integration
			double y_result[6];
//			Set initial condition (we only set this once, at the beginning of the integration
			y_result[0] = init_mag[0];
			y_result[1] = init_mag[1];
			y_result[2] = init_mag[2];
			y_result[3] = init_mag[3];
			y_result[4] = init_mag[4];
			y_result[5] = init_mag[5]; 
			#endif
			
			#ifdef SWEEP_FIELD_AND_CURRENT
			for(count = 0; count < 6; count++)
				y_result[count] = y_result_down[count];
			#endif
			
//			Do the actual integration!
//			Note that if INT_TIME_DEBUG is defined, we should save the starting time so that we can print out
//			debug messages if this integration step takes a while
			#ifdef INT_TIME_DEBUG
			time_t int_start_time;
			time(&int_start_time);
			int curr_step = 1;
			#endif
			
			#ifdef AVERAGE_MAGNETIZATION
			double running_weighted_sum = 0.0;
			double t_last = 0.0;
			#endif
			
			#ifdef TRACK_TSETTLING_X
			for(count = 0; count < 3; count++) {
				last_three_lyr1_max[count]=0;
				last_three_lyr1_pts[count]=0;
				last_three_lyr2_max[count]=0;
				last_three_lyr2_pts[count]=0;
			}
			
			tsettling_lyr1 = 0;
			tsettling_lyr2 = 0;
			#endif
			#ifdef TRACK_TSETTLING_Z
			tsettling_lyr1 = 0;
			tsettling_lyr2 = 0;
			#endif
			
			while(t < t_f) {
				//if INT_TIME_DEBUG defined, check that we've been integrating for less than INT_TIME_SECONDS seconds
				//if we haven't, print a debug message, and reset start_time
				#ifdef INT_TIME_DEBUG
				time_t int_curr_time;
				time(&int_curr_time);
				//printf("start was %d, curr is %d\n", (int)int_start_time, (int)int_curr_time);
				if(((int)int_curr_time) - ((int)int_start_time) > INT_TIME_SECONDS*curr_step) {
					fprintf(stderr, "DEBUG: Sweep down running on i = %d, k = %d, taken longer than %f seconds (t = %e, t_f = %e)\n", 
						i, k, INT_TIME_SECONDS * curr_step, t, t_f);
					curr_step++;
				}
				#endif
					
				//if we do a single-point solution, save the solution
				if((n==1) && (n_field==1)) {					 
					//If the number of evolutions is equal to the current return value allocated size, reallocate and ask for double the memory
					//I'm not sure why the -1 is, other than I'm scared of off-by-one errors
					if(num_evolutions == (retval_size-1)) {
						retval = realloc(retval, (retval_size * 2 * 4 * sizeof(double)));
						if(retval == NULL) {
							fprintf(stderr, "ERROR: Single-point solution output vector memory allocation failed, couldn't allocate %ld bytes. Maybe the machine ran out of memory?\n",
								(retval_size * 2 * 4 * sizeof(double)) );
							exit(-1);
						}
						retval_size *= 2;
					}
					//Save solution to retval
					retval[num_evolutions*4] = t;
					if(LOOK_AT_LAYER_2==1) {
						retval[num_evolutions*4 + 1] = y_result[3];
						retval[num_evolutions*4 + 2] = y_result[4];
						retval[num_evolutions*4 + 3] = y_result[5];
					}
					else {
						retval[num_evolutions*4 + 1] = y_result[0];
						retval[num_evolutions*4 + 2] = y_result[1];
						retval[num_evolutions*4 + 3] = y_result[2];
					}
				}
					
				int status = gsl_odeiv_evolve_apply (e, c, s, sys_ptr, &t, t_f, &h, y_result);
     			num_evolutions++;
     			
				if (status != GSL_SUCCESS) {
					printf("Integrator broke for some reason - help!\n");
					break;
				}
				
				//Enforce a minimum step size, if defined
				#ifdef INT_MIN_STEP_SZ
				if(h < INT_MIN_STEP_SZ) h = INT_MIN_STEP_SZ;
				#endif
				
				//Print vector length if MIN_VEC_LEN_PRINT is greater than the current vector length
				#ifdef MIN_VEC_LEN_PRINT
				double veclp = sqrt( pow(y_result[0],2) + pow(y_result[1],2) + pow(y_result[2],2) );
				if(veclp < MIN_VEC_LEN_PRINT) {
					fprintf(stderr, "Vector length of %f less than minimum %f, at t=%e\n", veclp, MIN_VEC_LEN_PRINT, t);
					//printf("%f\t%e\n", veclp, t);
				}
				#endif
				
				//Enforce a minimum vector length, if defined
				#ifdef INT_RENORM_VECTOR
				double vecl = sqrt( pow(y_result[0],2) + pow(y_result[1],2) + pow(y_result[2],2) );
				double vecl2 = sqrt( pow(y_result[3],2) + pow(y_result[4],2) + pow(y_result[5],2) );
				if(vecl < INT_RENORM_VECTOR_TOL || vecl > (1.0/INT_RENORM_VECTOR_TOL)) {
					y_result[0] = y_result[0] * (1.0 / vecl);
					y_result[1] = y_result[1] * (1.0 / vecl);
					y_result[2] = y_result[2] * (1.0 / vecl);
				}
				if(vecl2 < INT_RENORM_VECTOR_TOL || vecl2 > (1.0/INT_RENORM_VECTOR_TOL)) {
					y_result[3] = y_result[3] * (1.0 / vecl2);
					y_result[4] = y_result[4] * (1.0 / vecl2);
					y_result[5] = y_result[5] * (1.0 / vecl2);
				}
				#endif
				
				//If we're calculating GMR, calculate it
				#ifdef CALCULATE_GMR
				y_res_5 = (y_result[5] >= 0.999) ? 0.999 : ((y_result[5] <= -0.999) ? -0.999 : y_result[5]);
				y_res_2 = (y_result[2] >= 0.999) ? 0.999 : ((y_result[2] <= -0.999) ? -0.999 : y_result[2]);
				curr_value = -1*cos(acos( y_res_5 ) - acos( y_res_2 )); //calculate GMR
				#else
				curr_value = (LOOK_AT_LAYER_2==1) ? y_result[5] : y_result[2];
				#endif
				
				//if AVERAGE_MAGNETIZATION is defined, add current result to weighted sum
				#ifdef AVERAGE_MAGNETIZATION
				running_weighted_sum += (t-t_last)*curr_value;
				t_last = t;
				#endif
				
				#ifdef TRACK_TSETTLING_X
				//This is the code that looks and calculates the switching time
				last_three_lyr1_pts[0] = last_three_lyr1_pts[1];
				last_three_lyr1_pts[1] = last_three_lyr1_pts[2];
				last_three_lyr1_pts[2] = y_result[1];
				
				//Did we just find a maximum
				if(last_three_lyr1_pts[0] < last_three_lyr1_pts[1] && last_three_lyr1_pts[2] < last_three_lyr1_pts[1]) {
					last_three_lyr1_max[0] = last_three_lyr1_max[1];
					last_three_lyr1_max[1] = last_three_lyr1_max[2];
					last_three_lyr1_max[2] = y_result[1];
					 
					if(last_three_lyr1_max[1] >= TSETTLING_CRITERIA && last_three_lyr1_max[2] < TSETTLING_CRITERIA && tsettling_lyr1 == 0) {
						tsettling_lyr1 = t;
					}
				}

				//This is the code that looks and calculates the switching time
				//For the second layer (if the second layer is actually being calculated)
				if(numdims > 3) {
					last_three_lyr2_pts[0] = last_three_lyr2_pts[1];
					last_three_lyr2_pts[1] = last_three_lyr2_pts[2];
					last_three_lyr2_pts[2] = y_result[4];
				
					//Did we just find a maximum
					if(last_three_lyr2_pts[0] < last_three_lyr2_pts[1] && last_three_lyr2_pts[2] < last_three_lyr2_pts[1]) {
						last_three_lyr2_max[0] = last_three_lyr2_max[1];
						last_three_lyr2_max[1] = last_three_lyr2_max[2];
						last_three_lyr2_max[2] = y_result[4];
						 
						if(last_three_lyr2_max[1] >= TSETTLING_CRITERIA && last_three_lyr2_max[2] < TSETTLING_CRITERIA && tsettling_lyr2 == 0) {
							tsettling_lyr2 = t;
						}
					}
				}
				#endif
				
				//Below block gets executed if we are tracking TSETTLING based on a Z criterion
				#ifdef TRACK_TSETTLING_Z
				if(y_result[2] < TSETTLING_CRITERIA_Z && tsettling_lyr1==0) tsettling_lyr1=t;
				if(numdims > 3) {
					if(y_result[5] < TSETTLING_CRITERIA_Z && tsettling_lyr2==0) tsettling_lyr2=t;
				}
				#endif
			}
			
			if((n==1) && (n_field==1)) {
				printf("Did %d evolutions for single point solution\n", num_evolutions);
			}
			
			//Track the switching point
			#ifdef TRACK_SWITCHING_POINT
			if((i*n_field + k - 1) >= 0 && sd_init[i*n_field + k - 1] == 1) {
				if((prev_value_sd[i*n_field + k - 1] < 0.0 && curr_value > 0.0) || (prev_value_sd[i*n_field + k - 1] > 0.0 && curr_value < 0.0)) {
					sd_switch_hz = Hz[k];
					sd_switch_j = j[i];
				
					printf("Sweep down switch\t\t%e\t\t%e\n", sd_switch_hz, sd_switch_j);
					//printf("DEBUG: curr_value=%e, prev_value_sd=%e\n", curr_value, prev_value_sd);
				}
			}
			prev_value_sd[i*n_field + k] = curr_value;
			sd_init[i*n_field + k]=1;
			#endif
			
			//Print TSETTLING if either TSETTLING_X or TSETTLING_Z is defined for this run
			#ifdef TRACK_TSETTLING_X
			if(tsettling_lyr1 != 0)
				fprintf(stderr,"ALPHA: Found L1 (tsettling, hz)=\t%e\t%e\n", tsettling_lyr1, Hz[k]);
			if(tsettling_lyr2 != 0)
				fprintf(stderr,"ALPHA: Found L2 (tsettling, hz)=\t%e\t%e\n", tsettling_lyr2, Hz[k]);
			#endif
			#ifdef TRACK_TSETTLING_Z
			if(tsettling_lyr1 != 0)
				fprintf(stderr,"ALPHA: Found L1 (tsettling, hz)=\t%e\t%e\n", tsettling_lyr1, Hz[k]);
			if(tsettling_lyr2 != 0)
				fprintf(stderr,"ALPHA: Found L2 (tsettling, hz)=\t%e\t%e\n", tsettling_lyr2, Hz[k]);
			#endif
			
//			Store result in magsweepdown
			//magsweepdown[i*n_field + k] = y_result[2];
			
			#if (AVERAGE_MAGNETIZATION)
				magsweepdown[i*n_field+k] = running_weighted_sum; 
			#else
				magsweepdown[i*n_field + k] = curr_value;
			#endif
			
			#ifdef SWEEP_FIELD_AND_CURRENT
			for(count = 0; count < 6; count++)
				y_result_down[count] = y_result[count];
			#endif
			
//			Free the things we defined a while ago
			gsl_odeiv_evolve_free (e);
			gsl_odeiv_control_free (c);
			gsl_odeiv_step_free (s);
			free(param_list);
		} //end magsweepdown for loop
		
		#ifdef SWEEP_CURRENT
//		Make space for the output of the integration

//		Set initial condition (multiply one of them by -1, depending on what the "free" layer we're looking at is
		#if LOOK_AT_LAYER_2 == 0
		y_result[0] = init_mag[0];
		y_result[1] = init_mag[1];
		y_result[2] = (-1.0) * init_mag[2];
		y_result[3] = init_mag[3];
		y_result[4] = init_mag[4];
			#ifdef SWITCH_BOTH_LAYERS
				y_result[5] = (-1.0) * init_mag[5];
			#endif
			#ifndef SWITCH_BOTH_LAYERS
				y_result[5] = init_mag[5];
			#endif
		#elif LOOK_AT_LAYER_2 == 1
		y_result[0] = init_mag[0];
		y_result[1] = init_mag[1];
			#ifdef SWITCH_BOTH_LAYERS
				y_result[2] = (-1.0) * init_mag[2];
			#endif
			#ifndef SWITCH_BOTH_LAYERS
				y_result[2] = init_mag[2];
			#endif
		y_result[3] = init_mag[3];
		y_result[4] = init_mag[4];
		y_result[5] = (-1.0) * init_mag[5];
		#endif
		
		#endif //end restart from Y values
		
//		sweep up
		for(k = 0; k < n_field; k++) {
//			put together a parameter list...
			double ** param_list = (double **) malloc(27 * sizeof(double *));
			param_list[0] = init_mag;
			param_list[1] = Hz;
			param_list[2] = j;
			param_list[3] = &Beta_div_P_and_j;
			param_list[4] = Alpha;
			param_list[5] = Alpha_bar;
			param_list[6] = H_K1;
			param_list[7] = H_K2;
			param_list[8] = &thetaoffset;
			param_list[9] = (double *) ND;
			param_list[10] = (double *) NR;
			param_list[11] = (double *) NR_down;
			param_list[12] = Ms;
			param_list[13] = Gamma_bar;
			param_list[14] = (double *) &k; //yes, k is an int, so what
			param_list[15] = (double *) &i; //yes, i is an int, this may not be portable but should be fine for x86(-64)
			param_list[16] = (double*) &j_timedep_funs;
			param_list[17] = (double*) &hz_timedep_funs;
			param_list[18] = (double*) &temperature;
			param_list[19] = &dx;
			param_list[20] = &dy;
			param_list[21] = dz;
			param_list[22] = layer_position_z;
			param_list[23] = af_coupling_js;
			param_list[24] = pol;
			param_list[25] = &t_relax;
			param_list[26] = (double*)&sim_field_in_plane;
			
//			set up the ode system (most of this code is lifted from the GSL documentation and just re-written)
//			http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html
//			we don't give it a jacobian because we're going to use runge-kutta algorithm which doesn't need it
			
			gsl_odeiv_system* sys_ptr;
			
			//2D simulation case
			if(symtype == 6 || symtype == 9) {
				gsl_odeiv_system sys = {llg1, NULL, 3, param_list};
				numdims = 3;
				sys_ptr = &sys;
			}
			else if(symtype == 10) {
				gsl_odeiv_system sys = {llg1_6, NULL, 6, param_list};
				numdims = 6;
				sys_ptr = &sys;
			}
			//3D, FxFrFr or FrFrFx
			else if(symtype == 41 || symtype == 26) {
				gsl_odeiv_system sys = {llg1_fxfrfr, NULL, 6, param_list};
				numdims = 6;
				sys_ptr = &sys;
			}
			//If we're doing a simulation with the fixed layer on the BOTTOM reverse the direction of interactions
			//by reversing the down->up and up->down mutual demag tensors
			if(symtype == 9 || symtype == 41) {
				double* temp = param_list[10];
				param_list[10] = param_list[11];
				param_list[11] = param_list[10];
			}
			
//			define step type
//			So far, gsl_odeiv_step_rk4 seems to be the best, but gsl_odeiv_step_rk8pd isn't bad either. R-K (4,5) appears to suck.
			const gsl_odeiv_step_type * T = INTEGRATOR_TO_USE;
//			define the "stepper", give it the step type and # of dimensions
			gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, numdims);
//			Define controller
//			Ivan: Stephanie, please examine this - do we need to use the _y or _yp version? Also parameters to this controller, I just guessed.
//			see http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html
//			gsl_odeiv_control * c = gsl_odeiv_control_y_new (absolute_precision, relative_precision);
			gsl_odeiv_control * c = gsl_odeiv_control_standard_new (absolute_precision, relative_precision, 1.0, 1.0);
//			Allocate the function we're going to use for evolution
			gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (numdims);
//			Define the starting time for our integration
			double t = t_0;
//			Define the initial step size
			double h = t_f / INIT_STEP_SZ;
			
			#ifdef RESTART_Y_FROM_INIT_VALS
//			Make space for the output of the integration
			double y_result[6];
//			Set initial condition (multiply one of them by -1, depending on what the "free" layer we're looking at is
			#if LOOK_AT_LAYER_2 == 0
			y_result[0] = init_mag[0];
			y_result[1] = init_mag[1];
			y_result[2] = (-1.0) * init_mag[2];
			y_result[3] = init_mag[3];
			y_result[4] = init_mag[4];
				#ifdef SWITCH_BOTH_LAYERS
					y_result[5] = (-1.0) * init_mag[5];
				#endif
				#ifndef SWITCH_BOTH_LAYERS
					y_result[5] = init_mag[5];
				#endif
			#elif LOOK_AT_LAYER_2 == 1
			y_result[0] = init_mag[0];
			y_result[1] = init_mag[1];
				#ifdef SWITCH_BOTH_LAYERS
					y_result[2] = (-1.0) * init_mag[2];
				#endif
				#ifndef SWITCH_BOTH_LAYERS
					y_result[2] = init_mag[2];
				#endif
			y_result[3] = init_mag[3];
			y_result[4] = init_mag[4];
			y_result[5] = (-1.0) * init_mag[5];
			#endif
			
			#endif //end restart from Y values
			
			#ifdef SWEEP_FIELD_AND_CURRENT
			for(count = 0; count < 6; count++)
				y_result[count] = y_result_up[count];
			#endif
			
			//INT_TIME_DEBUG stuff to check that we don't run too long on a step
			#ifdef INT_TIME_DEBUG
			time_t int_start_time;
			time(&int_start_time);
			int curr_step = 1;
			#endif
			
			#ifdef AVERAGE_MAGNETIZATION
			double running_weighted_sum = 0.0;
			double t_last = 0.0;
			#endif
			
//			Do the actual integration!
			while(t < t_f) {
			
				//if INT_TIME_DEBUG defined, check that we've been integrating for less than INT_TIME_SECONDS seconds
				//if we haven't, print a debug message, and reset start_time
				#ifdef INT_TIME_DEBUG
				time_t int_curr_time;
				time(&int_curr_time);
				//printf("start was %d, curr is %d\n", (int)int_start_time, (int)int_curr_time);
				if(((int)int_curr_time) - ((int)int_start_time) > INT_TIME_SECONDS*curr_step) {
					fprintf(stderr, "DEBUG: Sweep up running on i = %d, k = %d, taken longer than %f seconds (t = %e, t_f = %e)\n", 
						i, k, INT_TIME_SECONDS * curr_step, t, t_f);
					curr_step++;
				}
				#endif
				
				int status = gsl_odeiv_evolve_apply (e, c, s, sys_ptr, &t, t_f, &h, y_result);
     
				if (status != GSL_SUCCESS) {
					printf("Integrator broke for some reason - help!\n");
					break;
				}
				
				//Enforce a minimum step size, if defined
				#ifdef INT_MIN_STEP_SZ
				if(h < INT_MIN_STEP_SZ) h = INT_MIN_STEP_SZ;
				#endif
				
				//Print vector length if MIN_VEC_LEN_PRINT is greater than the current vector length
				#ifdef MIN_VEC_LEN_PRINT
				double veclp = sqrt( pow(y_result[0],2) + pow(y_result[1],2) + pow(y_result[2],2) );
				if(veclp < MIN_VEC_LEN_PRINT) {
					fprintf(stderr, "Vector length of %f less than minimum %f, at t=%e\n", veclp, MIN_VEC_LEN_PRINT, t);
				}
				#endif
				
				#ifdef INT_RENORM_VECTOR
				double vecl = sqrt( pow(y_result[0],2) + pow(y_result[1],2) + pow(y_result[2],2) );
				double vecl2 = sqrt( pow(y_result[3],2) + pow(y_result[4],2) + pow(y_result[5],2) );
				if(vecl < INT_RENORM_VECTOR_TOL || vecl > (1.0/INT_RENORM_VECTOR_TOL)) {
					y_result[0] = y_result[0] * (1.0 / vecl);
					y_result[1] = y_result[1] * (1.0 / vecl);
					y_result[2] = y_result[2] * (1.0 / vecl);
				}
				if(vecl2 < INT_RENORM_VECTOR_TOL || vecl2 > (1.0/INT_RENORM_VECTOR_TOL)) {
					y_result[3] = y_result[3] * (1.0 / vecl2);
					y_result[4] = y_result[4] * (1.0 / vecl2);
					y_result[5] = y_result[5] * (1.0 / vecl2);
				}
				#endif
				
				//If we're calculating GMR, calculate it
				#ifdef CALCULATE_GMR
				y_res_5 = (y_result[5] >= 0.999) ? 0.999 : ((y_result[5] <= -0.999) ? -0.999 : y_result[5]);
				y_res_2 = (y_result[2] >= 0.999) ? 0.999 : ((y_result[2] <= -0.999) ? -0.999 : y_result[2]);
				curr_value = -1 * cos(acos( y_res_5 ) - acos( y_res_2 )); //calculate GMR
				#else
				curr_value = (LOOK_AT_LAYER_2==1) ? y_result[5] : y_result[2];
				#endif
				
				//if AVERAGE_MAGNETIZATION is defined, add current result to weighted sum
				#ifdef AVERAGE_MAGNETIZATION
				running_weighted_sum += (t-t_last)*curr_value;
				t_last = t;
				#endif
			}
			
			//printf("DEBUG: Last results for y[6] were {%.3e, %.3e, %.3e, %.3e, %.3e, %.3e}\n", y_result[0], y_result[1], y_result[2], y_result[3], y_result[4], y_result[5]);
			
			//fprintf(stderr, "Finished integrator!\n");
			
			//Track the switching point
			#ifdef TRACK_SWITCHING_POINT
			if((i*n_field + k - 1)>=0 && su_init[i*n_field + k - 1] == 1) {
				if((prev_value_su[i*n_field + k - 1] < 0.0 && curr_value > 0.0) || (prev_value_su[i*n_field + k - 1] > 0.0 && curr_value < 0.0)) {
					su_switch_hz = Hz[k];
					su_switch_j = j[i];
				
					printf("Sweep up switch\t\t%e\t\t%e\n", su_switch_hz, su_switch_j);
					//printf("DEBUG: curr_value=%e, prev_value_su=%e\n", curr_value, prev_value_su);
				}
			}
			su_init[i*n_field + k]=1;
			prev_value_su[i*n_field + k] = curr_value;
			#endif
			
//			Store result in magsweepup
			#if (AVERAGE_MAGNETIZATION)
				magsweepup[i*n_field+k] = running_weighted_sum; 
			#else
				magsweepup[i*n_field + k] = curr_value;
			#endif
			
			#ifdef SWEEP_FIELD_AND_CURRENT
			for(count = 0; count < 6; count++)
				y_result_up[count] = y_result[count];
			#endif

//			Free the things we defined a while ago
			gsl_odeiv_evolve_free (e);
			gsl_odeiv_control_free (c);
			gsl_odeiv_step_free (s);
			free(param_list);
		}//end magsweepup for loop				
	}

	time(&stop);
	fprintf(stderr, "About to free stuff!!!\n\n");
//	printf("Took %d seconds to run\n", (int) stop-start);

//	free memory stuff
	for(count = 0; count < 4; count++) {
		if(NR[count] != NULL) gsl_matrix_free(NR[count]);
		if(ND[count] != NULL) gsl_matrix_free(ND[count]);
		if(NR_down[count] != NULL) gsl_matrix_free(NR_down[count]);
	}
	
	for(count = 4; count < 6; count++) {
		if(NR[count] != NULL) gsl_matrix_free(NR[count]);
		if(NR_down[count] != NULL) gsl_matrix_free(NR_down[count]);
	}
	
	//If we're doing a single point solution, mark the end of the array with t = -1.0
	if((n==1) && (n_field==1)) {
		retval[4*num_evolutions] = -1.0;
	}
	
	fprintf(stderr, "Done freeing stuff.\n\n");
	
	#ifdef TRACK_SWITCHING_POINT
	printf("Switching point for sweep-down was\n\t%e along a field sweep\n\t%e along current\n\n", sd_switch_hz, sd_switch_j);
	printf("Switching point for sweep-up was\n\t%e along a field sweep\n\t%e along current\n\n", su_switch_hz, su_switch_j);
	printf("Center of switching window was \n\t%e along a field sweep\n\t%e along current\n\n", (su_switch_hz+sd_switch_hz)/2.0, (sd_switch_j+su_switch_j)/2.0);
	#endif
	
	return(retval);
}

/*This function wraps Hsweep to unpack parameter values for the multi-threaded implementation
param_ptr is a pointer to an array of double pointers, which contains all of the parameter values
packed inside of main
*/
void* Hsweep_wrapper(void* param_ptr) {
//	Re-cast the parameter as a (double**) value
	double** params = (double**) param_ptr;
		
//	Call Hsweep, note that we cleverly arrange the parameters of Hsweep to be in the same order as they were packed. Looking at this code it's a bit confusing, but if you look
//	at the code that calls Hsweep_wrapper() down in main() it makes a lot more sense what we are doing here
	Hsweep(*params[0], *params[1], params[2], params[3], params[4], params[5], *(int*)params[6], *params[7], *params[8], *(int*)params[9], *params[10], 
	*params[11], *(int*)params[12], *params[13], *(int*)params[14], *(int*)params[15], *params[16], *params[17], *params[18], params[19], params[20],
	params[21], params[22], params[23], params[24], *(int*)params[25], *(int*)params[26], params[27], params[28], *params[29]);
	
//	Doesn't return anything useful
	return(NULL);
}

/*List of parameters to pass in
(1) Theta offset (offset, in degrees, that the applied H-field is from the Z axis of the sample
(2) Percent that K2 (2nd order anisotropy) is of K1
(3) Alpha (gilbert's damping constant) for layer 1
(4-7) Layer1 through layer 4 saturation magnetization
(8-11) Layer 1 through layer 4 anisotropy (K1)
(12-15) Distance of layer 1 through layer 4 from bottom of sample, i.e. z=0
(16) Simulation type, lower 8 bits used. Each 2-bit pair encodes layers 1-4 (4 - LS bits) types. 00 means 
	disabled, 01 means fixed, 10 means free
(17-19) Minimum current density, maximum current density, number of current density points to plot against (units are A/m**2)
(20-22) Minimum applied field, maximum applied field, number of applied field points to plot against (units are teslas)
(23) Integration time, in seconds
(24) J time dependent functions, bit-wise OR'd (see documentation)
(25) H time dependent functions, bit-wise OR'd (see documentation)
(26-27) Sample size
(28) Temperature of the sample, for "thermal noise" - use 0 to ignore temperature-dependent noise
(29-40) Initial (x,y,z) magnetization for each of layers 1-4.
(41-44) Layer thickness for each of layers 1-4
(45-48) Polaization for each layer 
(49) Sample shape (integer, CURRENTLY UNUSED)
(50) Number of threads to use
(51-56) The af coupling for layer 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
(57-59) Alpha for layers 2,3,4
(60) Relaxation time before H and J get applied (seconds)
*/
int main(int argc, char ** argv) {
	//Check number of parameters on the command line, if not, print help message.
	if(argc != NUMARGS + 1) {
		printf("HELP:\nspinsim_gsl should be given %d parameters:\n", NUMARGS);
		printf("(1) Theta offset (offset, in degrees, that the applied H-field is from the Z axis of the sample\n(2) Percent that K2 (2nd order anisotropy) is of K1\
		\n(3) Alpha (gilbert's damping constant)\
			\n(4-7) Layer1 through layer 4 saturation magnetization\n(8-11) Layer 1 through layer 4 anisotropy (K1)\n\
			(12-15) Distance of layer 1 through layer 4 from bottom of sample, i.e. z=0\
			\n(16) Simulation type, lower 8 bits used. Each 2-bit pair encodes layers 1-4 (4 - LS bits) types. 00 means disabled, 01 means fixed, 10 means free\
			\n(17-19) Minimum current density, maximum current density, number of current density points to plot against (units are A/m**2)\
			\n(20-22) Minimum applied field, maximum applied field, number of applied field points to plot against (units are teslas)\n(23) Integration time, in seconds\
			\n(24) J time dependent functions, bit-wise OR'd (see documentation)\n(25) H time dependent functions, bit-wise OR'd (see documentation)\n\
			(26-27) Sample size dx, dy\
			\n(28) Temperature of the sample, for 'thermal noise' - use 0 to ignore temperature-dependent noise\n\
			(29-40) initial (x,y,z) magnetization for layers 1-4 respectively\
			(41-44) Layer thicknesses for layers 1-4\n(45-48) Polarization for each layer (CURRENTLY UNUSED)\n(49) Sample Shape (CURRENTLY UNUSED)\
			\n(50) Number of threads to use\n(51-56) The af coupling for layer 1-2, 1-3, 1-4, 2-3, 2-4, 3-4\n(57-59) Alpha for L2 - L4\n(60) Relaxation time before H/J application");
		exit(0);
	}
	
	int count;
	
	//unpack and parse in all of the command-line arguments, see comment at function heading
	double theta_offset = atof(argv[1]); //unpack theta_offset...etc...
	double per_K2 = atof(argv[2]);
	double Alpha[4];
	Alpha[0] = atof(argv[3]);
	//Unpack all 4 saturation magnetization values
	double Ms[4];
	for(count = 0; count < 4; count++) Ms[count] = atof(argv[4+count]);
	//Unpack all 4 anisotropy values
	double K1[4];
	for(count = 0; count < 4; count++) K1[count] = atof(argv[8+count]);
	//Unpack all 4 distances from z = 0
	double layer_position_z[4];
	for(count = 0; count < 4; count++) layer_position_z[count] = atof(argv[12+count]);
	//Unpack simulation type
	int sim_type = atoi(argv[16]);
	//Unpack Min/max/number of current density
	double j_min = atof(argv[17]);
	double j_max = atof(argv[18]);
	int n = atoi(argv[19]);
	//Unpack Min/max/number of applied field
	double hz_min = atof(argv[20]);
	double hz_max = atof(argv[21]);
	int n_field = atoi(argv[22]);
	
	double t_f = atof(argv[23]);
	unsigned int j_timedep_funs = atoi(argv[24]);
	unsigned int hz_timedep_funs = atoi(argv[25]);
	//Unpack sample size
	double dx = atof(argv[26]);
	double dy = atof(argv[27]);
	//Unpack Temperature
	double temperature = atof(argv[28]);
	//Unpack initial magnetizations
	double init_mag[12];
	for(count = 0; count < 12; count++) init_mag[count] = atof(argv[29+count]);
	//Unpack layer thicknesses
	double dz[4];
	for(count = 0; count < 4; count++) dz[count] = atof(argv[41+count]);
	//Unpack polarizations (NOT IMPLEMENT YET)
	double pol[4];
	for(count = 0; count < 4; count++) pol[count] = atof(argv[45+count]);
	//Unpack sample shape (NOT IMPLEMENT YET)
	;
	//read in # of threads to use.
	int numthreads = atoi(argv[50]);
	//Read in AF coupling constants
	double af_coupling_js[6];
	for(count = 0; count < 6; count++) af_coupling_js[count] = atof(argv[51+count]);
	Alpha[1]=atof(argv[57]);
	Alpha[2]=atof(argv[58]);
	Alpha[3]=atof(argv[59]);
	double t_relax;
	t_relax=atof(argv[60]);
	//[/end unpacking of parameters]
	
	//Declare and initialize where we store results from Hsweep()'s calculations
	//for these values we'll only use the result from the 0'th threads calculations
	double * Hz = (double*) malloc(sizeof(double) * n_field);
	double * j = (double*) malloc(sizeof(double) * n);
	//Save pointers of "garbage variables" to be free'd later, just so we can pass a "garbage" pointer to the 1st and later threads
	double** to_free_hz = (double**) malloc(sizeof(double*) * numthreads);
	double** to_free_j = (double**) malloc(sizeof(double*) * numthreads);
	//Note on threading: all threads will share this memory space. The threads are configured so they will
	//NOT crash into each other's data
	double * magsweepup = (double*) malloc(sizeof(double) * n * n_field);
	double * magsweepdown = (double*) malloc(sizeof(double) * n * n_field);
	//This pointer gets used if we do a single point (n==1, n_field==1) call ONLY
	//It stores a pointer to the returned memory space holding the x/y/z/t axis data
	double * retpointer_single;
	
	//make all of the magsweepup sections - we're going to pass these to the individual threads, and then re-combine them
	double ** magsweepup_threads = (double**)malloc(sizeof(double*) * numthreads);
	double ** magsweepdown_threads = (double**)malloc(sizeof(double*) * numthreads);
	for(count = 0; count < numthreads; count++) {
		magsweepup_threads[count] = (double*) malloc(sizeof(double) * n * n_field);
		magsweepdown_threads[count] = (double*) malloc(sizeof(double) * n * n_field);
	}
	
	//Keep track of the threads we've created
	int * thread_numbers = (int*) malloc(sizeof(int) * numthreads);
	//benchmarking
	time_t start, stop;
	//loop variables
	int count2;
	
	//Seed random number generator
	int tim;
	time((time_t*) &tim);
	srand(tim);
	
	//Load in time-dependent field/current function parameters from timedep function definition files, if necessary
	if(j_timedep_funs > 1 || hz_timedep_funs > 1) {
		if(LoadTDParams() != 0) {
			printf("ERROR: LoadTDParams() returned an error. Exiting...");
			exit(-1);
		}
	}
	
	//Threading stuff!
	//Allocate pointer to all pthread_t 
	pthread_t *threads = (pthread_t *) malloc(numthreads * sizeof(*threads));
	//Attributes for creating threads
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	//Specify that all threads are joinable
	pthread_attr_setdetachstate(&pthread_custom_attr,PTHREAD_CREATE_JOINABLE);
	
	time(&start);
	printf("Spinsim Hsweep Started!\n");

//	Single-threaded call
//	Note that if n and n_field are both 1 (single-point solution) multi-threading doesn't make sense so we single thread that anyway
	if( numthreads == 1 || ((n == 1) && (n_field == 1)) ) {
		retpointer_single = Hsweep(theta_offset, per_K2, Alpha, Ms, K1, 
		layer_position_z, sim_type, j_min, j_max, n, hz_min, 
		hz_max, n_field, t_f, j_timedep_funs, hz_timedep_funs, dx, dy, 
		temperature, init_mag, dz, j, Hz, magsweepup, magsweepdown, numthreads, 
		0, af_coupling_js, pol, t_relax);
	}
	
//	Multi-threaded call, means packing and passing parameters!
	else {
//		Launch all of numthreads threads
		for(count = 0; count < numthreads; count++) {
//			Allocate space and pack in all of the parameters that will be passed to the count-th thread
			double** params = (double**) malloc(30 * sizeof(double *));
//			Start packing the parameters!						
			params[0] = &theta_offset;
			params[1] = &per_K2;
			params[2] = Alpha;
			params[3] = Ms;
			params[4] = K1;
			params[5] = layer_position_z;
			params[6] = (double*) &sim_type;
			params[7] = &j_min;
			params[8] = &j_max;
			params[9] = (double*) &n;
			params[10] = &hz_min;
			params[11] = &hz_max;
			params[12] = (double*) &n_field;
			params[13] = &t_f;
			params[14] = (double*) &j_timedep_funs;
			params[15] = (double*) &hz_timedep_funs;
			params[16] = &dx;
			params[17] = &dy;
			params[18] = &temperature;
			params[19] = init_mag;
			params[20] = dz;
			params[25] = (double*) &numthreads;
			params[27] = af_coupling_js;
			params[28] = pol;
			params[29] = &t_relax;
						
			//If we're calling the first threads, pass it the "real" j and Hz. Otherwise, pass it the garbage variables
			//that will be ignored when the calculation is done
			if(count == 0) {
				params[22] = Hz;
				params[21] = j;
			}
			else {
				params[22] = (double*) malloc(n_field * sizeof(double));
				params[21] = (double*) malloc(n * sizeof(double));
				to_free_hz[count-1] = params[22];
				to_free_j[count-1] = params[21];
			}
			params[23] = magsweepup_threads[count];
			params[24] = magsweepdown_threads[count];
			
			//save thread number, and how many threads there are
			thread_numbers[count] = count;
			//Tell the thread what thread number it is
			params[26] = (double*) &thread_numbers[count];
			
			//call as an actual thread procedure
			pthread_create(&threads[count], &pthread_custom_attr, Hsweep_wrapper, (void*) params);
		}
	
		//Wait for all threads to finish (join)
		for (count=0; count < numthreads; count++)
		{
			pthread_join(threads[count],NULL);
		}

		//free garbage variables
		for(count = 0; count < numthreads-1; count++) {
			free(to_free_hz[count]);
			free(to_free_j[count]);
		}
	}
	
//	Combine all of the magsweepup results
//	Note that threads are given interleaved external current density lines to compute;

	//fprintf(stderr, "GOT HERE 1\n");

	if(numthreads > 1) {
		int h_count, j_count;
		 for(j_count = 0; j_count < n; j_count++){
			for(h_count = 0; h_count < n_field; h_count++) {
				//fprintf(stderr, "About to access j_count=%d, h_count=%d, addr=%d thread=%d (max=%d)\n", j_count, h_count, j_count*n_field + h_count, j_count % numthreads, n*n_field);
				magsweepup[j_count*n_field + h_count] = magsweepup_threads[(j_count%numthreads)][j_count*n_field + h_count];
				magsweepdown[j_count*n_field + h_count] = magsweepdown_threads[(j_count%numthreads)][j_count*n_field + h_count];
			}
		}
	}

	//fprintf(stderr, "GOT HERE 2\n");
	
	time(&stop);
	printf("100 percent Done!");
	printf("Total time elapsed was %d seconds\n", (int) (stop-start));
	
	//output files for Hsweep output
	//Put together all of the file names based on the parameters given
	FILE* output_fileup, *output_filedown, *j_output, *hz_output;
	char filename[128];
	
// If we are doing a single-point solution, dump...differently...
	if((n==1) && (n_field==1)) {
		//File where we will output (newline delimited) x-axis data
		sprintf(filename, "magsingle_x");
	#if LOOK_AT_LAYER_2 == 1
		strcat(filename, "2");
	#endif
		output_fileup = fopen(filename,"w");
		//File where we will output (newline delimited) y-axis data
		sprintf(filename, "magsingle_y");
	#if LOOK_AT_LAYER_2 == 1
		strcat(filename, "2");
	#endif
		output_filedown = fopen(filename,"w");
		//File where we will output (newline delimited) z-axis data
		sprintf(filename, "magsingle_z");
	#if LOOK_AT_LAYER_2 == 1
		strcat(filename, "2");
	#endif
		j_output = fopen(filename,"w");
		//File where we will output (newline delimited) time data
		sprintf(filename, "magsingle_t");
	#if LOOK_AT_LAYER_2 == 1
		strcat(filename, "2");
	#endif
		hz_output = fopen(filename,"w");
		
		//we'll reset this to 0 when we hit the end of the simulation
		int continues = 1;
		int curr_pos = 0;
		
		while(continues == 1) {
			//if we hit a time stamp of -1, this marks the end of our four-some list
			if(retpointer_single[curr_pos] == -1.0) {
				continues = 0;
			}
			//Put the current four-some to each file
			else {
				fprintf(hz_output, "%e\n", retpointer_single[curr_pos]);
				fprintf(output_fileup, "%e\n", retpointer_single[curr_pos+1]);
				fprintf(output_filedown, "%e\n", retpointer_single[curr_pos+2]);
				fprintf(j_output, "%e\n", retpointer_single[curr_pos+3]);
				curr_pos += 4;
			}
		}
	}
//	Otherwise, we're doing a multi-point solution. Only print the solution at t_f for this case
	else {	
		//We're going to put the magsweepup/down/j/Hz data in the following files
		sprintf(filename, "magsweepup_output");
		output_fileup = fopen(filename,"w");
	
		sprintf(filename, "magsweepdown_output");
		output_filedown = fopen(filename,"w");
	
		sprintf(filename, "magsweep_output_j");
		j_output = fopen(filename,"w");
	
		sprintf(filename, "magsweep_output_hz");
		hz_output = fopen(filename,"w");
	
		//Write out magsweepdown data
		for(count = 0; count < n; count++) {
			for(count2 = 0; count2 < n_field; count2++) {
				fprintf(output_filedown, "%.3e  ", magsweepdown[count*n_field + count2]);
			
				if(magsweepdown[count*n_field + count2] > 0.0) {
					fprintf(output_filedown, " ");
				}
			}
			fprintf(output_filedown, "\n");
		}
	
		//write out magsweepup data
		for(count = 0; count < n; count++) {
			for(count2 = 0; count2 < n_field; count2++) {
				fprintf(output_fileup, "%.3e  ", magsweepup[count*n_field + count2]);
			
				if(magsweepup[count*n_field + count2] > 0.0) {
					fprintf(output_fileup, " ");
				}
			}
			fprintf(output_fileup, "\n");
		}
	
		//write out j and hz
		for(count = 0; count < n; count++) {
			fprintf(j_output, "%.3e\n", j[count]);
		}
		for(count = 0; count < n_field; count++) {
			fprintf(hz_output, "%.3e\n", Hz[count]);
		}
	}
	
//	Close up all of our files
	fclose(output_fileup);
	fclose(output_filedown);
	fclose(j_output);
	fclose(hz_output);
	
	return(0);
}

