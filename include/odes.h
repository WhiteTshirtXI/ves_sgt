/* ORDINARY DIFFERENTIAL EQUATIONS
 *  Governing equations dx(t)/dt = f(x,t) for the boundary-value problem.
 *
 * REFERENCES
 *	See supporting documentation for the derivation of the equations.
 */

#ifndef ODES_H
#define ODES_H

/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "sys.h"

/* PROTOTYPES */
void odesS0(double a, double b, double t, double *p, double *x, double *dx);

/* IMPLEMENTATIONS */

void odes(char &id, int order, double a, double b, double t, 
          double *p, double *x, double *dx){
	
	// calculate differential
	if      (order == 0){
		if      (id == 'S' ) odesS0(a, b, t, p, x, dx);
		//else if (id == "SB")
		//else if (id == "LR")
		//else if (id == "LS")
		//else if (id == "LR")
	}
	else if (order == 1){
	}
}

// Short vesicle, zeroth order (without bending elasticity)
void odesS0(double a, double b, double t, double *p, double *x, double *dx){
	// Assume memory has been allocated for dx.

	int M = 14;
	double D = b - a;

	// variables
	double p0r = x[ 0];	// rear pressure
	double s0r = x[ 1];	// rear slope
	double k0r = x[ 2];	// rear modified height
	double p0f = x[ 3];	// front pressure
	double s0f = x[ 4];	// front slope
	double k0f = x[ 5];	// front modified height
	double T0  = x[ 6];	// membrane tension
	double Q0  = x[ 7];	// leakback flux per unit circumference
	double C0r = x[ 8];	// rear intercept
	double C0f = x[ 9];	// front intercept
	double f0r = x[10];	// rear wall force density
	double f0f = x[11];	// front wall force density
	double F0r = x[12];	// rear wall force
	double F0f = x[13];	// front wall force

	double h0r = k0r + 0.5*s0r*s0r;	// rear height
	double h0f = k0f + 0.5*s0f*s0f;	// front height

	// differentials
	double dp0r = -(-6.0/h0r/h0r + 12.0*Q0/h0r/h0r/h0r);
	double ds0r = -(-1.0 - p0r/T0);
	double dk0r = -(s0r*(1.0 + ds0r));
	double dp0f =   -6.0/h0f/h0f + 12.0*Q0/h0f/h0f/h0f ;
	double ds0f =   -1.0 - p0f/T0 ;
	double dk0f =   s0f*(1.0 + ds0f) ;
	double dT0  = 0.0;
	double dQ0  = 0.0;
	double dC0r = 0.0;
	double dC0f = 0.0;
	double df0r = -(4.0/h0r - 6.0*Q0/h0r/h0r);
	double df0f =   4.0/h0f - 6.0*Q0/h0f/h0f ;
	double dF0r = 0.0;
	double dF0f = 0.0;

	dx[ 0] = D*dp0r;
	dx[ 1] = D*ds0r;
	dx[ 2] = D*dk0r;
	dx[ 3] = D*dp0f;
	dx[ 4] = D*ds0f;
	dx[ 5] = D*dk0f;
	dx[ 6] = D*dT0 ;
	dx[ 7] = D*dQ0 ;
	dx[ 8] = D*dC0r;
	dx[ 9] = D*dC0f;
	dx[10] = D*df0r;
	dx[11] = D*df0f;
	dx[12] = D*dF0r;
	dx[13] = D*dF0f;
}

// Short vesicle, first order (without bending elasticity)
void odesS1(){
}

// Short vesicle, zeroth order (with bending elasticity)
void odesSB0(){
}

// Short vesicle, first order (with bending elasticity)
void odesSB1(){
}

// Long vesicle (rear), zeroth order (without bending elasticity)
void odesLR0(){
}

// Long vesicle (front), zeroth order (without bending elasticity)
void odesLF0(){
}

// Long vesicle (rear), first order (without bending elasticity)
void odesLR1(){
}

// Long vesicle (front), first order (without bending elasticity)
void odesLF1(){
}



#endif
