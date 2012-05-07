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

//Todo now: 
//Plot!
//GUI
//Clean up and document

//Defines for constants
#define T_FIN 2.5e-8 	//time to end integration
#define EM_ES 650000	//free layer saturation magnetization
#define EM_AR_ES 600000	//pinned layer saturation magnetization
#define HZ_MAX 1e6		//Maximum (and minimum) applied external H-field
#define J_MAX 1e11		//Maximum (and minimum) applied external current density

//Defines for general program parameters
#define NOISINESS 5		//How talkative? 1 is "as loud as possible", a large number makes it quiet (like 10000)
//#define ELLIPSE		//undefine to use ellipse instead of circular model

//Defines for integrator
#define ABS_PREC 1e-6			//Absolute precision for integrator to aim for
#define REL_PREC 1e-6			//Relative precision for integrator to aim for
#define INIT_STEP_SZ 100000.0	//Define the proportion that the initial step will be of t_f, i.e. it will be t_f/INIT_STEP_SZ

//Defines for the testing function (main)
#define ENN_FIELD 30	//# of conditions to plot H-field against
#define ENN 30			//# of conditions to plot current density against
#define PERC_K2 0.5		//Percentage that K1 is of K2
#define THETA_OFF 0.05	//Theta offset, in degrees

//Defines for multi-threading
#define NUM_THREADS 1	//define # of threads to use

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

/* llg1() is a function that is the "differential equation" for the ODE solver
It takes parameters as defined by GSL odeiv system (http://www.gnu.org/software/gsl/manual/html_node/Defining-the-ODE-System.html)
t = the time position to use
y[] = the values of y at which our equation is being computed
dydt[] = the array in which to store the output
params = a pointer to a memory space of 19 pointers that point to the following variables:
	{p, Hz, j, Beta_div_P_and_j, Alpha, Alpha_bar, H_K1, H_K2, thetaoffset, ND, NR, Ms, MRs, m_x, m_y, m_z, Gamma_bar, k, i}
The function will return the values of dy/dt for each of the three dimensions STORED INTO dydt[]
The actual return value is GSL_SUCCESS for a valid execution */
int llg1(double t, const double y[], double dydt[], void * params) {
//	First, extract all values from params parameter
//	Ivan: We make a poor non-portable assumption that pointers to all types have same length
//	But it's probably good enough for most architectures
	double ** param_list = (double **) params;
	double * p = param_list[0];
	double * Hz = param_list[1];
	double * j = param_list[2];
	double Beta_div_P_and_j = *param_list[3];
	double Alpha = *param_list[4];
	double Alpha_bar = *param_list[5];
	double H_K1 = *param_list[6];
	double H_K2 = *param_list[7];
	double thetaoffset = *param_list[8];
	gsl_matrix* ND = (gsl_matrix*) param_list[9];
	gsl_matrix* NR = (gsl_matrix*) param_list[10];
	double Ms = *param_list[11];
	double MRs = *param_list[12];
	double mx_0 = *param_list[13];
	double my_0 = *param_list[14];
	double mz_0 = *param_list[15];
	double Gamma_bar = *param_list[16];
	int k = *((int*)param_list[17]);
	int i = *((int*)param_list[18]);
	double HP_K1 = *param_list[19];

//	Now we'll define and calculate variable's we'll actually use
	double Beta = Beta_div_P_and_j*polarisation(y[0],y[1],y[2],p[3],p[4],p[5])*j[i];
	double Beta_bar1 = Alpha*Beta/( 1 + pow(Alpha,2) );
	double Beta_bar2 = Beta/(1+pow(Alpha,2));

	double Betap = -Beta_div_P_and_j*polarisation(y[3],y[4],y[5],p[0],p[1],p[2])*j[i];
        double Beta_bar1p = Alpha*Betap/( 1 + pow(Alpha,2) );
        double Beta_bar2p = Betap/(1+pow(Alpha,2));

//	Re-define y as y_mat, a matrix with the values as different columns
	gsl_matrix* y_mat = gsl_matrix_alloc(3,1);
	gsl_matrix_set(y_mat, 0, 0, y[0]);
	gsl_matrix_set(y_mat, 1, 0, y[1]);
	gsl_matrix_set(y_mat, 2, 0, y[2]);
//	Re-define p as p_mat, a matrix with the values as different columns
	gsl_matrix* p_mat = gsl_matrix_alloc(3,1);
	gsl_matrix_set(p_mat, 0, 0, y[3]);
	gsl_matrix_set(p_mat, 1, 0, y[4]);
	gsl_matrix_set(p_mat, 2, 0, y[5]);

//      Free Layer
//	Define H, the applied field, as a gsl_matrix, remember that it's (rows, columns)
	gsl_matrix* H = gsl_matrix_alloc(3,1);
	gsl_matrix_set(H, 0, 0, 0.0);
	gsl_matrix_set(H, 1, 0, 0.0);
	gsl_matrix_set(H, 2, 0, Hz[k]);
//	Anisotropic field
	gsl_matrix* Hani = gsl_matrix_alloc(3,1);
	gsl_matrix_set(Hani, 0, 0, H_K1*y[0]*sin(thetaoffset));
	gsl_matrix_set(Hani, 1, 0, 0.0);
	gsl_matrix_set(Hani, 2, 0, H_K1*y[2]*cos(thetaoffset)+H_K2*y[2]*( pow(y[1],2) + pow(y[2],2) ) );
//	Self-demag field "Hani"
	gsl_matrix* HD = gsl_matrix_alloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND, y_mat, 0, HD);
	gsl_matrix_scale(HD, (-1.0 * Ms));
//	Mutual demag field "HR"
	gsl_matrix* HR = gsl_matrix_alloc(3,1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR, p_mat, 0, HR);
	gsl_matrix_scale(HR, (-1.0 * MRs));
//	Effective field is the sum of the anisotropy field (Hani), the self (HD)
//	and mutual (HR) demag. fields, and the external field (H)
	gsl_matrix* Heff = gsl_matrix_calloc(3,1); //calloc sets initial values to 0
	gsl_matrix_add(Heff, H);
	gsl_matrix_add(Heff, Hani);
	gsl_matrix_add(Heff, HR);
	gsl_matrix_add(Heff, HD);

//      Pinned Layer
	gsl_matrix* HPani = gsl_matrix_alloc(3,1);
        gsl_matrix_set(Hani, 0, 0, 0.0);
        gsl_matrix_set(Hani, 1, 0, 0.0);
        gsl_matrix_set(Hani, 2, 0, HP_K1*y[5]);
	gsl_matrix* HPD = gsl_matrix_alloc(3,1);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ND, p_mat, 0, HD);
        gsl_matrix_scale(HD, (-1.0 * MRs));
	gsl_matrix* HPR = gsl_matrix_alloc(3,1);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, NR, y_mat, 0, HR);
        gsl_matrix_scale(HR, (-1.0 * Ms));
        gsl_matrix* HPeff = gsl_matrix_calloc(3,1); //calloc sets initial values to 0                                                                                                                   
        gsl_matrix_add(HPeff, H);
        gsl_matrix_add(HPeff, HPani);
        gsl_matrix_add(HPeff, HPR);
        gsl_matrix_add(HPeff, HPD);

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
	dydt[0] = -1.0 * Gamma_bar * (y[1]*gsl_matrix_get(Heff, 2, 0)-y[2]*gsl_matrix_get(Heff, 1, 0)) - Alpha_bar*((y[0]*y[0]-1)*gsl_matrix_get(Heff, 0, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 1, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 2, 0)) + Beta_bar1 * (y[1]*y[5] - y[2]*y[4]) - Beta_bar2 * ((y[0]*y[0]-1)*y[3]+y[0]*y[1]*y[4]+y[0]*y[2]*y[5]);
	dydt[1] = -1.0 * Gamma_bar*(y[2]*gsl_matrix_get(Heff, 0, 0)-y[0]*gsl_matrix_get(Heff, 2, 0))-Alpha_bar*((y[1]*y[1]-1)*gsl_matrix_get(Heff, 1, 0)+y[0]*y[1]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 2, 0))+Beta_bar1*(y[2]*y[3]-y[0]*y[5]) - Beta_bar2*((y[1]*y[1]-1)*y[4]+y[0]*y[1]*y[3]+y[1]*y[2]*y[5]);
	dydt[2] = -1.0 * Gamma_bar*(y[0]*gsl_matrix_get(Heff, 1, 0)-y[1]*gsl_matrix_get(Heff, 0, 0))-Alpha_bar*((y[2]*y[2]-1)*gsl_matrix_get(Heff, 2, 0)+y[0]*y[2]*gsl_matrix_get(Heff, 0, 0)+y[1]*y[2]*gsl_matrix_get(Heff, 1, 0))+Beta_bar1*(y[0]*y[4]-y[1]*y[3])-Beta_bar2*((y[2]*y[2]-1)*y[5]+y[0]*y[2]*y[3]+y[1]*y[2]*y[4]);
	
	// Ref Layer Dynamics
	dydt[3] = -1.0 * Gamma_bar*(y[4]*gsl_matrix_get(HPeff, 2, 0)-y[5]*gsl_matrix_get(HPeff, 1, 0)) - Alpha_bar*((y[3]*y[3]-1)*gsl_matrix_get(HPeff, 0, 0)+y[3]*y[4]*gsl_matrix_get(HPeff, 1, 0)+y[3]*y[5]*gsl_matrix_get(HPeff, 2, 0))+Beta_bar1p*(y[4]*y[2]-y[5]*y[1])-Beta_bar2p*((y[3]*y[3]-1)*y[0]+y[3]*y[4]*y[1]+y[3]*y[5]*y[2]);
	dydt[4] = -1.0*Gamma_bar*(y[5]*gsl_matrix_get(HPeff, 0, 0)-y[3]*gsl_matrix_get(HPeff, 2, 0))-Alpha_bar*((y[4]*y[4]-1)*gsl_matrix_get(HPeff, 1, 0)+y[3]*y[4]*gsl_matrix_get(HPeff, 0, 0)+y[4]*y[5]*gsl_matrix_get(HPeff, 2, 0))+Beta_bar1p*(y[5]*y[0]-y[3]*y[2])-Beta_bar2p*((y[4]*y[4]-1)*y[1]+y[3]*y[4]*y[0]+y[4]*y[5]*y[2]);
	dydt[5] = -1.0*Gamma_bar*(y[3]*gsl_matrix_get(HPeff, 1, 0)-y[4]*gsl_matrix_get(HPeff, 0, 0))-Alpha_bar*((y[5]*y[5]-1)*gsl_matrix_get(HPeff, 2, 0)+y[3]*y[5]*gsl_matrix_get(HPeff, 0, 0)+y[4]*y[5]*gsl_matrix_get(HPeff, 1, 0))+Beta_bar1p*(y[3]*y[1]-y[4]*y[0])-Beta_bar2p*((y[5]*y[5]-1)*y[2]+y[3]*y[5]*y[0]+y[4]*y[5]*y[1]);
	
	//Checked that this appears to behave correctly 2009-01-03
	
	//free things we've used, mostly all of the matrices
	gsl_matrix_free(Heff);
	gsl_matrix_free(HR);	
	gsl_matrix_free(Hani);
	gsl_matrix_free(HD);
	gsl_matrix_free(H);
	gsl_matrix_free(p_mat);
	gsl_matrix_free(y_mat);

	return(GSL_SUCCESS);
}

/* Hsweep() calculates the magnetization of a free layer across a 2D plane of external field & applied current conditions
Parameters:
n = the number of current values to evaluate at
nfield = the number of applied h-field values to evaluate at
percK2 = the percentage that K2 is of K1
thetaoff = the offset of the applied h-field vs the z-axis of the structure
Hz = a 1D-array with size nfield containing the h-field conditions applied to the structure
magsweepup = a 2D array of size (n, nfield) for the magnetization on upsweeping calculation
magsweepdown = a 2D array of size (n, nfield) for the magnetization on downsweeping calculation
j = a 1D-array with size n of the applied current conditions
thread_no = the thread number this is. 0 = first thread
*/
int Hsweep(int n, int nfield, double percK2, double thetaoff, double * Hz, double * magsweepup, double * magsweepdown, double * j, int thread_no) {
//	save start and stop times to benchmark the program speed
	time_t start, stop;
	time(&start);

//	Initial conditions in time and magnetization
	double mz_0 = .9; 
	double mx_0 = 0; 
	double my_0 = .1;
//	Magnitude of the magnetization vector
	double m = sqrt( pow((mx_0),2) + pow((my_0),2) + pow((mz_0),2) );
//	Make the magnetization vector unit sized
	mx_0 = mx_0 / m; 
	my_0 = my_0 / m;
	mz_0 = mz_0 / m; 

//	Integration times
	double t_0 = 0;
	double t_f = T_FIN;

//	Mutual Demagnetization Tensor
	double X = 0;
	double Y = 0;
	double Z = 7;
//	circle
	double dx=50;
	double dy=50;
	double  dz=3;
//	ellipse
	#ifdef ELLIPSE
		dx = 100; dy = 50; dz = 3;
	#endif
	
	gsl_matrix * NR = gsl_matrix_alloc(3,3);
	demagtensor(X,Y,Z,dx,dy,dz,NR);
//	Note: For infinite lateral size use NR_xx = 0; NR_yy = 0; NR_zz = 0;

//	Self Demagnetization Matrix
//	X,Y,Z,dx,dy,dz were already declared, this is minus style points
	X = 0; Y = 0; Z = 0;
//	circle
	dx=50; dy=50; dz=3;
//	ellipse
//	dx = 100; dy = 50; dz = 3;
	gsl_matrix * ND = gsl_matrix_alloc(3,3);
	demagtensor(X,Y,Z,dx,dy,dz,ND);
//	Note: For infinite lateral size use ND_xx = 0; ND_yy = 0; ND_zz = 1;

//	Define the normal LLG parameters for the problem
	double Gamma = 2.21e5; 			//Gyromagnetic ratio, 2.2127616 x 1e5 s^(-1).A^(-1).m
	double Alpha = 0.01; 			//Gilbert's damping factor, dimensionless
	double d = 3e-9;               //Free layer thickness, m
	double Ms = EM_ES; 			//Free layer saturation magnetization, A/m
	double K1 = 3e5;               //Free layer anisotropy constant, J/m^3
	double K2 = (percK2/100)*K1;   //Free layer anisotropy constant
	double mu0 = 4*PI*1e-7; 		//Magnetic permeability of vaccuum, H/m
	double e = 1.602e-19; 			//Electron charge, -1.602176 x 1e-19 C
	double hbar = 1.054e-34; 		//Reduced Planck constant, 1.054571 x 1e-34 J.s

//      Define params for the reference layer - for reasonable dynamics, make 1.5x K1
	double Kp1 = 3e5*10;          // pinned layer anisotropy

//	Other constants
	double MRs = EM_AR_ES;                //Reference layer saturation magnetization, A/m
	double px = 0; 
	double py = 0; 
	double pz = 1;		 //Pinned layer initial state
	double p[3] = {0, 0, 1};
		
//	Slonczewski's spin transfert torque coefficient Beta=P*j*Beta_div_P_and_j
	double Beta_div_P_and_j = (Gamma*hbar)/(2*d*mu0*Ms*e);

//	Some useful parameters
	double H_K1 = 2*K1/mu0/Ms;                   //Free layer anisotropy field, A/m
	double H_K2 = 4*K2/mu0/Ms;                   //Free layer anisotropy field, A/m
	double HP_K1 = 2*Kp1/mu0/Ms;                        //Pinned layer anisotropy field, A/m
	double norm = H_K1 / sqrt( pow(H_K1,2)+pow(H_K2,2) );
	H_K1 = norm*H_K1;
	H_K2 = norm*H_K2;
	//H_K2 = 0; //This doesn't feel right, but it was originally uncommented
	double Gamma_bar = Gamma/(1+pow(Alpha,2) );
	double Alpha_bar = Alpha*Gamma/( 1 + pow(Alpha,2) );

//	Initialize our applied field vector
//	Applied field will go from (-HZ_MAX) to (HZ_MAX) with nfield steps in between
//	TODO: This would probably be better if the space did not have to center about 0
	make_linspace(HZ_MAX,nfield,Hz);

//	Current Density Vector
//	Applied current density will go from (-J_MAX) to (J_MAX) with nfield steps in between
	make_linspace(J_MAX,n,j); 

//	Offset for easy axis of free FM layer
	double thetaoffset = thetaoff*PI/180;
	
//	Loop variables
	int i, k;
	
	for(i = thread_no; i < n; i+=NUM_THREADS) {
	  //if((i-thread_no) % (NOISINESS*NUM_THREADS) == 0) {
			printf("Currently on i = %d (of %d) iteration (thread %d)\n", i, n, thread_no);
			//}

//		sweep down
		for(k = (nfield-1); k >= 0; k--) {
//			put together a parameter list...
			double ** param_list = (double **) malloc(19 * sizeof(double *));
			param_list[0] = p;
			param_list[1] = Hz;
			param_list[2] = j;
			param_list[3] = &Beta_div_P_and_j;
			param_list[4] = &Alpha;
			param_list[5] = &Alpha_bar;
			param_list[6] = &H_K1;
			param_list[7] = &H_K2;
			param_list[8] = &thetaoffset;
			param_list[9] = (double *) ND;
			param_list[10] = (double *) NR;
			param_list[11] = &Ms;
			param_list[12] = &MRs;
			param_list[13] = &mx_0;
			param_list[14] = &my_0;
			param_list[15] = &mz_0;
			param_list[16] = &Gamma_bar;
			param_list[17] = (double *) &k; //yes, k is an int, so what
			param_list[18] = (double *) &i; //yes, i is an int, this may not be portable but should be fine for x86(-64)
			param_list[19] = &HP_K1;

//			set up the ode system (most of this code is lifted from the GSL documentation and just re-written)
//			http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html
//			we don't give it a jacobian because we're going to use runge-kutta algorithm which doesn't need it
			gsl_odeiv_system sys = {llg1, NULL, 6, param_list};
//			define step type
//			So far, gsl_odeiv_step_rk4 seems to be the best, but gsl_odeiv_step_rk8pd, gsl_odeiv_step_rk2 isn't bad either. R-K (4,5) appears to suck.
			const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
//			define the "stepper", give it the step type and # of dimensions
			gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 6);
//			Define controller
//			Ivan: Stephanie, please examine this - do we need to use the _y or _yp version? Also parameters to this controller, I just guessed.
//			see http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html
			gsl_odeiv_control * c = gsl_odeiv_control_y_new (ABS_PREC, REL_PREC);
//			Allocate the function we're going to use for evolution
			gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (6);
//			Define the starting time for our integration
			double t = t_0;
//			Define the initial step size
			double h = t_f / INIT_STEP_SZ;
//			Make space for the output of the integration
			double y_result[6];
//			Set initial condition
			y_result[0] = mx_0;
			y_result[1] = my_0;
			y_result[2] = mz_0;
			y_result[3] = px;
                        y_result[4] = py;
                        y_result[5] = pz;
			
//			Do the actual integration!
			while(t < t_f) {
				int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t_f, &h, y_result);
			
     
				if (status != GSL_SUCCESS) {
					printf("Integrator broke for some reason - help!\n");
					break;
				}
			}
			
//			Store result in magsweepdown
			magsweepdown[i*nfield + k] = y_result[2];

//			Free the things we defined a while ago
			gsl_odeiv_evolve_free (e);
			gsl_odeiv_control_free (c);
			gsl_odeiv_step_free (s);
			free(param_list);
		} //end magsweepdown for loop
		
//		sweep up
		for(k = 0; k < nfield; k++) {
//			put together a parameter list...
			double ** param_list = (double **) malloc(19 * sizeof(double *));
			param_list[0] = p;
			param_list[1] = Hz;
			param_list[2] = j;
			param_list[3] = &Beta_div_P_and_j;
			param_list[4] = &Alpha;
			param_list[5] = &Alpha_bar;
			param_list[6] = &H_K1;
			param_list[7] = &H_K2;
			param_list[8] = &thetaoffset;
			param_list[9] = (double *) ND;
			param_list[10] = (double *) NR;
			param_list[11] = &Ms;
			param_list[12] = &MRs;
			param_list[13] = &mx_0;
			param_list[14] = &my_0;
			param_list[15] = &mz_0;
			param_list[16] = &Gamma_bar;
			param_list[17] = (double *) &k; //yes, k is an int, so what
			param_list[18] = (double *) &i; //yes, i is an int, this may not be portable but should be fine for x86(-64)
			param_list[19] = &HP_K1;

//			set up the ode system (most of this code is lifted from the GSL documentation and just re-written)
//			http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html
//			we don't give it a jacobian because we're going to use runge-kutta algorithm which doesn't need it
			gsl_odeiv_system sys = {llg1, NULL, 6, param_list};
//			define step type
//			So far, gsl_odeiv_step_rk4 seems to be the best, but gsl_odeiv_step_rk8pd isn't bad either. R-K (4,5) appears to suck.
			const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;
//			define the "stepper", give it the step type and # of dimensions
			gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 6);
//			Define controller
//			Ivan: Stephanie, please examine this - do we need to use the _y or _yp version? Also parameters to this controller, I just guessed.
//			see http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html
			gsl_odeiv_control * c = gsl_odeiv_control_y_new (ABS_PREC, REL_PREC);
//			Allocate the function we're going to use for evolution
			gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (6);
//			Define the starting time for our integration
			double t = t_0;
//			Define the initial step size
			double h = t_f / INIT_STEP_SZ;
//			Make space for the output of the integration
			double y_result[6];
//			Set initial condition
			y_result[0] = mx_0;
			y_result[1] = my_0;
			y_result[2] = (-1.0) * mz_0;
			y_result[0] = px;
                        y_result[1] = py;
                        y_result[2] = pz;
			
//			Do the actual integration!
			while(t < t_f) {
				int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t_f, &h, y_result);
     
				if (status != GSL_SUCCESS) {
					printf("Integrator broke for some reason - help!\n");
					break;
				}
			}
			
			//fprintf(stderr, "Finished integrator!\n");
			
//			Store result in magsweepdown
			magsweepup[i*nfield + k] = y_result[2];

//			Free the things we defined a while ago
			gsl_odeiv_evolve_free (e);
			gsl_odeiv_control_free (c);
			gsl_odeiv_step_free (s);
			free(param_list);
		}//end magsweepup for loop				
	}

	time(&stop);
//	printf("Took %d seconds to run\n", (int) stop-start);
//	free memory stuff
	gsl_matrix_free(NR);
	gsl_matrix_free(ND);
}

/*This function wraps Hsweep to unpack parameter values for the multi-threaded implementation
*/
int Hsweep_wrapper(void* param_ptr) {
	double** params = (double**) param_ptr;
	int n = *(int*) params[0];
	int n_field = *(int*) params[1];
	double per_k2 = *(double*) params[2];
	double theta_offset = *(double*) params[3];
	double* hz_ptr = (double*) params[4];
	double* magsweepup = (double*) params[5];
	double* magsweepdown = (double*) params[6];
	double* j = (double*) params[7];
	int thread_number = *(int*) params[8];
	
	Hsweep(n, n_field, per_k2, theta_offset, hz_ptr, magsweepup, magsweepdown, j, thread_number);
}

int main() {
	//for these values we'll only use the result from the 0'th threads calculations
	double Hz[ENN_FIELD];
	double j[ENN];
	//Save pointers of "garbage variables" to be free'd later
	double* to_free_hz[NUM_THREADS];
	double* to_free_j[NUM_THREADS];
	//Note on threading: all threads will share this memory space. The threads are configured so they will
	//NOT crash into each other's data
	double magsweepup[ENN][ENN_FIELD];
	double magsweepdown[ENN][ENN_FIELD];
	//So that we can pass these as "parameters" by "reference"
	double per_k2 = PERC_K2;
	double theta_offset = THETA_OFF;
	int n = ENN;
	int n_field = ENN_FIELD;
	int thread_numbers[NUM_THREADS];
	//benchmarking
	time_t start, stop;
	
	int count;
	
	//Allocate pointer to all pthread_t 
	pthread_t *threads = (pthread_t *) malloc(NUM_THREADS * sizeof(*threads));
	//Attributes for creating threads
	pthread_attr_t pthread_custom_attr;
	pthread_attr_init(&pthread_custom_attr);
	//Specify that all threads are joinable
	pthread_attr_setdetachstate(&pthread_custom_attr,PTHREAD_CREATE_JOINABLE);
	
	time(&start);

//	Single-threaded call
	if(NUM_THREADS == 1) {
		Hsweep(ENN, ENN_FIELD, PERC_K2, THETA_OFF, Hz, (double*)magsweepup, (double*)magsweepdown, j, 0);
	}
	else {
		//Launch all of NUM_THREADS threads
		for(count = 0; count < NUM_THREADS; count++) {
			//Allocate space and pack in all of the parameters that will be passed to the count-th thread
			double** params = (double**) malloc(9 * sizeof(double *));
			params[0] = (double*) &n;
			params[1] = (double*) &n_field;
			params[2] = (double*) &per_k2;
			params[3] = (double*) &theta_offset;
			params[5] = (double*) magsweepup;
			params[6] = (double*) magsweepdown;
		
			//For j and Hz, we'll only save the results from the 0-th thread's calculations
			//as the arrays will be the same across all threads
			if(count == 0) {
				params[4] = (double*) Hz;
				params[7] = (double*) j;
			}
			else {
				params[4] = (double*) malloc(ENN_FIELD * sizeof(double));
				params[7] = (double*) malloc(ENN * sizeof(double));
				to_free_hz[count-1] = (double*)params[4];
				to_free_j[count-1] = (double*)params[7];
			}
		
			thread_numbers[count] = count;
			params[8] = (double*) &thread_numbers[count];
		
			pthread_create(&threads[count], &pthread_custom_attr, Hsweep_wrapper, (void*)	params);
		}
	
		//Wait for all threads to finish
		for (count=0; count < NUM_THREADS; count++)
		{
			pthread_join(threads[count],NULL);
		}

		//free garbage variables
		for(count = 0; count < NUM_THREADS-1; count++) {
			free(to_free_hz[count]);
			free(to_free_j[count]);
		}
	}
	
	time(&stop);
	printf("Total time elapsed was %d seconds\n", (int) (stop-start));
	
	int count2;
	
	//output files for Hsweep output
	//TODO: make these names depend on parameters
	FILE* output_fileup, *output_filedown, *j_output, *hz_output;
	output_fileup = fopen("magsweep_outputup_30_30_0-5_0-05","w");
	output_filedown = fopen("magsweep_outputdown_30_30_0-5_0-05","w");
	j_output = fopen("magsweep_output_30_30_j","w");
	hz_output = fopen("magsweep_output_30_30_hz","w");
	
	//Write out magsweepdown data
	printf("Printing magsweepdown:\n");
	//fprintf(output_file, "Printing magsweepdown:\n");
	for(count = 0; count < ENN; count++) {
		for(count2 = 0; count2 < ENN_FIELD; count2++) {
			printf("%.3e  ", magsweepdown[count][count2]);
			fprintf(output_filedown, "%.3e  ", magsweepdown[count][count2]);
			if(magsweepdown[count][count2] > 0.0) {
				printf(" ");
				fprintf(output_filedown, " ");
			}
		}
		printf("\n");
		fprintf(output_filedown, "\n");
	}
	
	printf("\n\n");
	//fprintf(output_file, "\n\n");
	
	//write out magsweepup data
	printf("Printing magsweepup:\n");
	//fprintf(output_file, "Printing magsweepup:\n");
	for(count = 0; count < ENN; count++) {
		for(count2 = 0; count2 < ENN_FIELD; count2++) {
			printf("%.3e  ", magsweepup[count][count2]);
			fprintf(output_fileup, "%.3e  ", magsweepup[count][count2]);
			if(magsweepup[count][count2] > 0.0) {
				printf(" ");
				fprintf(output_fileup, " ");
			}
		}
		printf("\n");
		fprintf(output_fileup, "\n");
	}
	
	//write out j and hz
	for(count = 0; count < ENN; count++) {
		fprintf(j_output, "%.3e\n", j[count]);
	}
	for(count = 0; count < ENN_FIELD; count++) {
		fprintf(hz_output, "%.3e\n", Hz[count]);
	}
	
//	Close up all of our files
	fclose(output_fileup);
	fclose(output_filedown);
	fclose(j_output);
	fclose(hz_output);
	
	return(0);
}
