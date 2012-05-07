/* 	spin_helpers.c
	This file provides various helper functions for the spinsim.c top-level file.
	Note: This was written and tested on Linux 2.6.18 x86-64 with GSL version 1.12, it is recommended to use
	the same versions to prevent compatibility issues.
	Note: If compiling using gcc, be sure to use the '-lgsl -lgslcblas -lm' options on the command line.
	Written I. Yulaev 2010-01-04 (iyulaev@ucsd.edu)
	Engineered S. Moyerman (smoyerman@ucsd.edu)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Includes from GSL library
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

//This puts in a correction for the demag tensor calculation; we divide dz by 2 to give us
//half-height for the elliptic cylinder. This produces much better agreement to the example
//cited in the appendix of Xiao, Zangwill, Stiles PRB 72,014446 (2005)
#define CORRECT_DEMAGTENSOR

//Use an artificially even (about 90-degrees) g(theta) function
//#define SYMMETRIC_G_THETA

//Use g(theta) = g(0) for all theta
//#define CONSTANT_G_THETA

//Use LLG "tunneling" form of g(theta)
//don't use with Xiao form
//#define TUNNELING_G_THETA

//Use alternate form of tunneling current, where eta(theta) = sin(theta)
//#define TUNNELING_G_THETA2

//For the Xiao form, B0, B1, Qp, and Qn are hard-coded in the functions polarisation() and
//polarisationpp(). 

//math.h should include PI but whatever
#ifndef PI
#define PI (3.141592653589793L)
#endif

/* Returns the Newell_f output given parameters a,b,c
TODO: Add some relevant description.
Translated from Matlab "Newell_f.m" 2009-12-31 */
double newell_f(double a, double b, double c) {
	if(a == 0) {
		a = pow(10.0,(-20.0));
	}
	if(b == 0) {
		b = pow(10.0,(-20.0));
	}
	if(c == 0) {
		c = pow(10.0,(-20.0));
	}
	
//	calculate Newellf
	double output;
	
	double a_2 = pow(a,2);
	double b_2 = pow(b,2);
	double c_2 = pow(c,2);
	
	output = 0.5*b*(c_2 - a_2) * asinh( b/sqrt(a_2+c_2) ) + 0.5 * c * (b_2-a_2) * asinh( c/sqrt(a_2+b_2) ) - a * b * c * atan((b*c)/(a*sqrt(a_2+b_2+c_2)))+(2*a_2-b_2-c_2)/6*sqrt( a_2+b_2+c_2 );

	return(output);
}

/* Calculates f2(), called from f1()
Used by demagtensor() to put together a tensor
Translated from Matlab "f2.m" 2010-01-01 */
double f2(double a, double b, double c) {
	double returned = newell_f(a,b,c)-newell_f(a,0,c)-newell_f(a,b,0)+newell_f(a,0,0);
	return(returned);
}


/* Calculates f1(), called by f(),
Used by demagtensor() to put together a tensor
Translated from Matlab "f1.m" 2010-01-01 */
double f1(double x, double y, double z, double dx, double dy, double dz) {
	double returned = f2(x,y,z)-f2(x,y-dy,z)-f2(x,y,z-dz)+f2(x,y-dy,z-dz);
	return(returned);
}

/* Calculates f(), called by Nii(),
Used by demagtensor() to put together a tensor
Translated from Matlab "f.m" 2010-01-01 */
double f(double x, double y, double z, double dx, double dy, double dz) {
	double returned = f1(x,y+dy,z+dz,dx,dy,dz)-f1(x,y,z+dz,dx,dy,dz)-f1(x,y+dy,z,dx,dy,dz)+f1(x,y,z,dx,dy,dz);
	return(returned);
}

/* Calculates nii(), called by demagtensor(),
Used by demagtensor() to put together a tensor
Translated from Matlab "Nii.m" 2010-01-01 */
double nii(double x, double y, double z, double dx, double dy, double dz) {
	double returned = (1/(4*PI*dx*dy*dz))*(2*f(x,y,z,dx,dy,dz)-f(x+dx,y,z,dx,dy,dz)-f(x-dx,y,z,dx,dy,dz));
	return(returned);
}

/* Returns the Newell_g output given parameters a,b,c
TODO: Add some relevant description.
Translated from Matlab "Newell_g.m" 2010-01-01 */
double newell_g(double a, double b, double c) {
	if(a == 0) {
		a = pow(10.0,(-20.0));
	}
	if(b == 0) {
		b = pow(10.0,(-20.0));
	}
	if(c == 0) {
		c = pow(10.0,(-20.0));
	}
	
//	calculate Newellg
	double output;
	
	double a_2 = pow(a,2);
	double b_2 = pow(b,2);
	double c_2 = pow(c,2);
	double c_3 = pow(c,3);
	
	output = a*b*c*asinh(c/sqrt(a_2+b_2))+b/6*(3*c_2-b_2)*asinh(a/sqrt(b_2+c_2))+a/6*(3*c_2-a_2)*asinh(b/sqrt(a_2+c_2))-c_3/6*atan((a*b)/(c*sqrt(a_2+b_2+c_2)))-c*b_2/2*atan((a*c)/(b*sqrt(a_2+b_2+c_2)))-c*a_2/2*atan((b*c)/(a*sqrt(a_2+b_2+c_2)))-a*b*sqrt(a_2+b_2+c_2)/3;

	return(output);
}

/* Calculates g2(), called from g1()
Used by demagtensor() to put together a tensor
Translated from Matlab "g2.m" 2010-01-01 */
double g2(double a, double b, double c) {
	double returned = newell_g(a,b,c)-newell_g(a,b,0.0);
	return(returned);
}


/* Calculates g1(), called by g(),
Used by demagtensor() to put together a tensor
Translated from Matlab "g1.m" 2010-01-01 */
double g1(double x, double y, double z, double dx, double dy, double dz) {
	double returned = g2(x+dx,y,z+dz)-g2(x+dx,y,z)-g2(x,y,z+dz)+g2(x,y,z);
	return(returned);
}

/* Calculates g(), called by Nij(),
Used by demagtensor() to put together a tensor
Translated from Matlab "g.m" 2010-01-01 */
double g(double x, double y, double z, double dx, double dy, double dz) {
	double returned = g1(x,y,z,dx,dy,dz)-g1(x,y-dy,z,dx,dy,dz)-g1(x,y,z-dz,dx,dy,dz)+g1(x,y-dy,z-dz,dx,dy,dz);
	return(returned);
}

/* Calculates nij(), called by demagtensor(),
Used by demagtensor() to put together a tensor (TODO: describe better)
Translated from Matlab "Nij.m" 2010-01-01 */
double nij(double x, double y, double z, double dx, double dy, double dz) {
	double returned = (1/(4*PI*dx*dy*dz))*(g(x,y,z,dx,dy,dz)-g(x-dx,y,z,dx,dy,dz)-g(x,y+dy,z,dx,dy,dz)+g(x-dx,y+dy,z,dx,dy,dz));
	return(returned);
}

/* This function puts together a demagnetizing tensor given the parameters
The return value is status (0 means good)
The matrix is put into the pointer '* output_tensor'
Translated from Matlab "demagtensor.m" 2010-01-01 */
int demagtensor(double x, double y, double z, double dx, double dy, double dz, gsl_matrix* output_tensor) {
//	The demagnetizing tensor is stored as a matrix of values that looks like:
//	( Nxx  Nxy  Nxz )
//	( Nyx  Nyy  Nyz )
//	( Nzx  Nzy  Nzz )
//	We use gsl_matrix_set to set the values in the matrix (already allocated for us) accordingly

	#ifdef CORRECT_DEMAGTENSOR
	dz /= 2;
	#endif

	gsl_matrix_set(output_tensor, 0, 0, nii(x,y,z,dx,dy,dz));
	gsl_matrix_set(output_tensor, 0, 1, nij(x,y,z,dx,dy,dz));
	gsl_matrix_set(output_tensor, 0, 2, nij(x,z,y,dx,dz,dy));
	gsl_matrix_set(output_tensor, 1, 0, nij(y,x,z,dy,dx,dz));
	gsl_matrix_set(output_tensor, 1, 1, nii(y,x,z,dy,dx,dz));
	gsl_matrix_set(output_tensor, 1, 2, nij(y,z,x,dy,dz,dx));
	gsl_matrix_set(output_tensor, 2, 0, nij(z,x,y,dz,dx,dy));
	gsl_matrix_set(output_tensor, 2, 1, nij(z,y,x,dz,dy,dx));
	gsl_matrix_set(output_tensor, 2, 2, nii(z,y,x,dz,dy,dx));

	return(0);
}

/*Polarisation computes the polarisation of something or other (TODO)
Called by the differential equation function
Note that there is one version for Xiao-based polarization calculation,
and another version for the Slonczewski form
Translated from Matlab "polarisation.m" 2010-01-02 */
#ifdef POLARISATION_XIAO
double polarisation(double m1_x, double m1_y, double m1_z, double m2_x, double m2_y, double m2_z) {
	double B0=1;
	double B1=0.5;
	double Qp=1;
	double Qn=-0.2; //-0.1;

	double cosTheta = m1_x*m2_x+m1_y*m2_y+m1_z*m2_z;
	
	#ifdef SYMMETRIC_G_THETA
	cosTheta = abs(cosTheta);
	#endif
	
	#ifdef CONSTANT_G_THETA
	double returned = Qp/(B0+B1*1)+Qn/(B0-B1*1);
	return(returned);
	#endif
	
	double returned = Qp/(B0+B1*cosTheta)+Qn/(B0-B1*cosTheta);
	return(returned);
}

//This function gets called to compute the eta(theta) function for computing spin torque transfer
//Aside from the magnetization of the two layers it also takes the polarization of the two layers.
//This IGNORES polarization and just does the same thing as polarisation() for POLARISATION_XIAO
double polarisationpp(double m1_x, double m1_y, double m1_z, double m2_x, double m2_y, double m2_z, double p1, double p2) {
	double B0=1;
	double B1=0.5;
	double Qp=1;
	double Qn=-0.2; //-0.1;

	double cosTheta = m1_x*m2_x+m1_y*m2_y+m1_z*m2_z;
	
	#ifdef SYMMETRIC_G_THETA
	cosTheta = abs(cosTheta);
	#endif
	
	#ifdef CONSTANT_G_THETA
	double returned = Qp/(B0+B1*1)+Qn/(B0-B1*1);
	return(returned);
	#endif
	
	double returned = Qp/(B0+B1*cosTheta)+Qn/(B0-B1*cosTheta);
	return(returned);
}
#endif

#ifdef POLARISATION_SLONCZEWSKI
//This function gets called when the parameter 'P' is NOT given; we use a default value of 35%
double polarisation(double m1_x, double m1_y, double m1_z, double m2_x, double m2_y, double m2_z) {
	double P = 0.7;

	double cosTheta = m1_x*m2_x+m1_y*m2_y+m1_z*m2_z;
	
	#ifdef SYMMETRIC_G_THETA
	cosTheta = abs(cosTheta);
	#endif
	
	#ifdef CONSTANT_G_THETA
	return(-4+(pow((1+P),3)*(3+1)/(4*pow(P,1.5))));
	#endif
	
	#ifdef TUNNELING_G_THETA
	return( P/(2*(1+P*P*cosTheta)) );
	#endif
	
	#ifdef TUNNELING_G_THETA2
	return(sin(acos(cosTheta)));
	#endif
	
	double returned =-4+(pow((1+P),3)*(3+cosTheta)/(4*pow(P,1.5)));
	return(1.0/returned);
}

//This function gets called when the parameter 'P' is given
double polarisationp(double m1_x, double m1_y, double m1_z, double m2_x, double m2_y, double m2_z, double P) {

	double cosTheta = m1_x*m2_x+m1_y*m2_y+m1_z*m2_z;
	
	#ifdef SYMMETRIC_G_THETA
	cosTheta = abs(cosTheta);
	#endif
	
	#ifdef CONSTANT_G_THETA
	return(-4+(pow((1+P),3)*(3+1)/(4*pow(P,1.5))));
	#endif
	
	#ifdef TUNNELING_G_THETA
	return( P/(2*(1+P*P*cosTheta)) );
	#endif
	
	#ifdef TUNNELING_G_THETA2
	return(sin(acos(cosTheta)));
	#endif
	
	//returned is the same as 'g' in the Slonczewski Polarization Model
	double returned =-4+(pow((1+P),3)*(3+cosTheta)/(4*pow(P,1.5)));
	return(1.0/returned);
}

//This function gets called to compute the eta(theta) function for computing spin torque transfer
//Aside from the magnetization of the two layers it also takes the polarization of the two layers.
double polarisationpp(double m1_x, double m1_y, double m1_z, double m2_x, double m2_y, double m2_z, double p1, double p2) {
	
	double cosTheta = m1_x*m2_x+m1_y*m2_y+m1_z*m2_z;
	
	#ifdef SYMMETRIC_G_THETA
	cosTheta = abs(cosTheta);
	#endif
	
	#ifdef CONSTANT_G_THETA
	return((4*p1*sqrt(p2))/((1+p1)*(1+p1) * (1+p2) * (3+1) - 16*p1*sqrt(p2)));
	#endif
	
	#ifdef TUNNELING_G_THETA
	return( p1/(2*(1+p1*p2*cosTheta)) );
	#endif
	
	#ifdef TUNNELING_G_THETA2
	return(sin(acos(cosTheta)));
	#endif
	
	return((4*p1*sqrt(p2))/((1+p1)*(1+p1) * (1+p2) * (3+cosTheta) - 16*p1*sqrt(p2)));
}

#endif

/*This function creates a double linspace and sticks it into output
The linspace will go from (-1 * magnitude) to magnitude with 'size' positions
TODO: Allow this to set the end points arbitrarily
Written Ivan Yulaev 2010-01-02
*/
void make_linspace(double min, double max, int size, double * output) {
	int count;
	
	double magnitude = (max - min)/2;
	
	if(size == 1) {
		output[0] = min;
	}
	else {
		for(count = 0; count < size; count++) {
			output[count] = min + (max-min)*(((double)count)/((double)(size-1)));
			//printf("Linspace putting %e for the %dth value, sum of %e and %e\n", output[count], count, min, (max-min)*(((double)count)/((double)size)-1.0));
		}
	}
}

/* This function simply calculates the cosine of the angle between two vectors
Parameters: x1, y1, z1 is the cartesian representation of the first vector.
x2, y2, z2 is the cartesian representation of the second vactor
Returns cosine of the angle between the two
Used in llg1 to calculate AF coupling */
double cos_angle_bw_two_vectors(double x1, double y1, double z1, double x2, double y2, double z2) {
	double one_abs = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
	double two_abs = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));
	double dot_p = (x1*x2) + (y1*y2) + (z1*z2);
	
	return(dot_p / (one_abs * two_abs));
}

/*main here for testing purposes only
int main() {
	int count;
	gsl_matrix* test_m = gsl_matrix_alloc(3,3);

	demagtensor(0.3, 0.5, 0.7, 0.03, 0.05, 7e-2, test_m);

	for(count = 0; count < 3; count++) {
		printf("%f\t%f\t%f\n", gsl_matrix_get(test_m, count, 0), gsl_matrix_get(test_m, count, 1), gsl_matrix_get(test_m, count, 2));
	}
	
		double * blah = NULL;
	
	printf("The size of a double pointer is %d\n", sizeof(blah));
	printf("The size of a int pointer is %d\n", sizeof(int *));
	printf("The size of a gsl_matrix pointer is %d\n", sizeof(gsl_matrix*));
	
	gsl_matrix* y_mat = gsl_matrix_alloc(3,1);
	gsl_matrix_set(y_mat, 0, 0, 2);
	gsl_matrix_set(y_mat, 1, 0, 4);
	gsl_matrix_set(y_mat, 2, 0, 6);
	
	gsl_matrix* result = gsl_matrix_alloc(3,1);
	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, test_m, y_mat, 0, result);
	
	printf("Result looks like:\n%f\n%f\n%f\n", gsl_matrix_get(result, 0, 0), gsl_matrix_get(result, 1, 0), gsl_matrix_get(result, 2, 0));

	return(0);
}*/
