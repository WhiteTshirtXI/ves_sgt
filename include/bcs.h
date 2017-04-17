/* BOUNDARY CONDITIONS
 *  Boundary conditions A*x(0) + B*x(1) = c for the boundary-value problem.
 *
 * REFERENCES
 *	See supporting documentation for the derivation of the equations.
 */

#ifndef BCS_H
#define BCS_H

/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <math.h>

/* PROTOTYPES */
void bcsS0(int N, int M, double a, double b, double *p, double *x, double *A, double *B, double *c);

/* IMPLEMENTATIONS */

void bcs(char id, int order, int N, int M, double a, double b, double *p, double *x,
         double *A, double *B, double *c){
	bcsS0(N, M, a, b, p, x, A, B, c);
}

// Short vesicle, zeroth order (without bending elasticity)
void bcsS0(int N, int M, double a, double b, double *p, double *x, double *A, double *B, double *c){
	// Assume memory has been allocated for A, B, c.

	double D = b - a;

	// variables
	double T0  = x[ 6];	// membrane tension
	double Q0  = x[ 7];	// leakback flux per unit circumference
	double C0r = x[ 8];	// rear intercept
	double C0f = x[ 9];	// front intercept
	double F0r = x[12];	// rear wall force
	double F0f = x[13];	// front wall force

	// parameters
	double d = p[0];		// midsection length
	double b1 = b - d;

	// boundary condition coefficients
	for (int i = 0; i < M; i++){
		for (int j = 0; i < M; i++){
			A[i*M+j] = 0.0;
			B[i*M+j] = 0.0;
		}
		c[i] = 0.0;
	}

	A[ 0*M+ 0] = 1.0;	A[ 0*M+ 3] = -1.0;		// midpoint continuity of pressure
	A[ 1*M+ 1] = 1.0;	A[ 1*M+ 4] = -1.0;		// midpoint continuity of slope
	A[ 2*M+ 2] = 1.0;	A[ 2*M+ 5] = -1.0;		// midpoint continuity of height
	A[ 3*M+ 2] = 1.0;	                  		// midpoint height
	A[ 4*M+10] = 1.0;	A[ 4*M+11] = -1.0;		// midpoint continuity of wall force density
	A[ 5*M+10] = 1.0;			   							// midpoint wall force
	
	B[ 6*M+ 0] = 1.0;	B[ 6*M+ 6] =  2.0;		// rear pressure BC
	B[ 7*M+ 1] = 1.0;	                  		// rear slope BC
	B[ 8*M+ 2] = 1.0;	B[ 8*M+ 8] = -1.0;		// rear height BC
	B[ 9*M+ 3] = 1.0;	B[ 9*M+ 6] =  2.0;		// front pressure BC
	B[10*M+ 4] = 1.0;	                  		// front slope BC
	B[11*M+ 5] = 1.0;	B[11*M+ 9] = -1.0;		// front height BC
	B[12*M+10] = 1.0;	B[12*M+12] = -1.0;		// rear wall force BC
	B[13*M+11] = 1.0;	B[13*M+13] = -1.0;		// front wall force BC

	c[ 3] = 1.0;													// midpoint height
	c[ 7] = -b1;													// rear slope BC
	c[10] =  b1;													// front slope BC

	// add truncation error terms to rear and front BCs
	double b2, b4, b6;
	double ep0, es0, eh0, ef0, ek0;
	double C0;
	for (int i = 0; i < 2; i++){
		if (i == 0) { b1 = -b + d; b2 = b1*b1; b4 = b2*b2; b6 = b2*b4; C0 = C0r;}
		if (i == 1) { b1 =  b - d; b2 = b1*b1; b4 = b2*b2; b6 = b2*b4; C0 = C0f;}

		// pressure
		ep0 =   8.0/(b1*b2) 
		    -  48.0*(2.0*Q0 + 2.0*C0)/(5.0*b1*b4) 
		    +  64.0/(T0*b6) 
		    + 288.0*C0*(2.0*Q0 + C0)/(7.0*b1*b6) 
		    - 768.0*(2.0*Q0 + 2.0*C0)/(5.0*T0*b2*b6);
		
		// slope
		es0 =   4.0/(T0*b2) 
		    -  12.0*(2.0*Q0 + 2.0*C0)/(5.0*T0*b2*b2) 
		    +  64.0/(5.0*T0*T0*b1*b4) 
		    +  48.0*C0*(2.0*Q0 + C0)/(7.0*T0*b2*b4) 
		    - 768.0*(2.0*Q0 + 2.0*C0)/(35.0*T0*T0*b1*b6) 
		    -  16.0*(-64.0 + 5.0*C0*C0*T0*T0*(6.0*Q0 + 2.0*C0))/(15.0*T0*T0*T0*b2*b6);
		
		// height
		eh0 = - 4.0/(T0*b1)
		    +   4.0*(2.0*Q0 + 2.0*C0)/(5.0*T0*b1*b2)
		    -  16.0/(5.0*T0*T0*b2*b2) 
		    -  48.0*C0*(2.0*Q0 + C0)/(35.0*T0*b1*b4) 
		    + 128.0*(2.0*Q0 + 2.0*C0)/(35.0*T0*T0*b2*b4) 
		    +  16.0*(-64.0 + 5.0*C0*C0*T0*T0*(6.0*Q0 + 2.0*C0))/(105.0*T0*T0*T0*b1*b6);
		
		// wall force density
		ef0 = - 8.0/(b1)
		    +   4.0*(6.0*Q0 + 4.0*C0)/(3.0*b1*b2)
		    -  16.0/(T0*b2*b2) 
		    -  16.0*C0*(6.0*Q0 + 2.0*C0)/(5.0*b1*b4);
		
		// modified height
		ek0 = eh0 - b1*es0 - 0.5*es0*es0;

		if (i == 0){
			c[ 6] += ep0;
			c[ 7] += es0;
			c[ 8] += ek0;
		}
		if (i == 1){
			c[ 9] += ep0;
			c[10] += es0;
			c[12] += ek0;
		}
	}
}

// Short vesicle, first order (without bending elasticity)
void bcsS1(){
}

// Short vesicle, zeroth order (with bending elasticity)
void bcsSB0(){
}

// Short vesicle, first order (with bending elasticity)
void bcsSB1(){
}

// Long vesicle (rear), zeroth order (without bending elasticity)
void bcsLR0(){
}

// Long vesicle (front), zeroth order (without bending elasticity)
void bcsLF0(){
}

// Long vesicle (rear), first order (without bending elasticity)
void bcsLR1(){
}

// Long vesicle (front), first order (without bending elasticity)
void bcsLF1(){
}



#endif
