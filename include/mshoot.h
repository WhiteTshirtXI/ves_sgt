/* MULTIPLE SHOOTING METHOD
 *  Solve the nonlinear system of ODEs dx(t)/dt = f(x,t) by the multiple shooting method.
 *
 * REFERENCES
 *  Stoer, J & Bulirsch, R. Introduction to Numerical Analysis (3rd ed).
 *		Springer. 2000. pp. 302-316, 559-561.
 *  Moin, P. Fundamentals of Engineering Numerical Analysis (2nd ed). 
 *		Cambridge University Press. 2010. pp. 64-70.
 */

#ifndef MSHOOT_H
#define MSHOOT_H

/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "odes.h"
#include "bcs.h"

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Given the guess s for the solution vector x(t;p), compute the function vector F to be zeroed 
 * and Jacobian matrix DF for Newton-Raphson iteration. 
 *   N : half the number of shooting segments for 2*N+1 shooting points  
 *   M : number of variables (order of the differential system)
 */
void func(int N, int M, double *t, double *p, double *s, double *F, double *DF){

	/*--- INITIALIZATION  ---*/
	double ds = 0.000001;

	// Allocate memory
	double *x     = (double*) calloc(        M, sizeof(double));
	double *dx    = (double*) calloc(        M, sizeof(double));

	double *u     = (double*) calloc( 2*N+1   , sizeof(double));	// unit domain
	double *u_r   = (double*) calloc( 2*N     , sizeof(double));	// from u = 0 to 0.5, in N intervals
	double *u_f   = (double*) calloc( 2*N     , sizeof(double));	// from u = 1 to 0.5, in N intervals
	
	double *sig   = (double*) calloc(        M, sizeof(double));	// test (input) vector
	double *sig0  = (double*) calloc(        M, sizeof(double));
	double *sigM  = (double*) calloc(        M, sizeof(double));
	double *sigp  = (double*) calloc(        M, sizeof(double));
	double *s_r   = (double*) calloc(   N   *M, sizeof(double));
	double *s_f   = (double*) calloc(   N   *M, sizeof(double));
	double *s_m   = (double*) calloc(        M, sizeof(double));
	
	double *chi   = (double*) calloc(        M, sizeof(double));	// output vector
	double *chip  = (double*) calloc(        M, sizeof(double));
	double *x_r   = (double*) calloc(   N   *M, sizeof(double));
	double *x_f   = (double*) calloc(   N   *M, sizeof(double));

	double *r     = (double*) calloc(        M, sizeof(double));	// boundary condition vector

	double *F     = (double*) calloc((2*N+1)*M, sizeof(double));	// vector to be zeroed
	double *F_r   = (double*) calloc((  N-1)*M, sizeof(double));
	double *F_f   = (double*) calloc((  N-1)*M, sizeof(double));
	double *F_m   = (double*) calloc(     2 *M, sizeof(double));
	double *F_rc  = (double*) calloc(        M, sizeof(double));

	double *DF    = (double*) calloc((2*N+1)*M*(2*N+1)*M, sizeof(double));	// Jacobian matrix
	double *DF_r  = (double*) calloc(   N   *M*(2*N+1)*M, sizeof(double));
	double *DF_f  = (double*) calloc(   N   *M*(2*N+1)*M, sizeof(double));
	double *DF_rc = (double*) calloc(        M*(2*N+1)*M, sizeof(double));

	double *I     = (double*) calloc(        M        *M, sizeof(double)); // identity matrix
	double *G     = (double*) calloc(        M        *M, sizeof(double)); // partial derivatives
	double *A     = (double*) calloc(        M        *M, sizeof(double)); // boundary conditions
	double *B     = (double*) calloc(        M        *M, sizeof(double));
	double *c     = (double*) calloc(                  M, sizeof(double));
	
	// Domain decomposition
	double a = t[    0];	// lower boundary
	double b = t[2*N+1];	// upper boundary

	for (int i = 0; i < 2*N+1; i++) 
		u[i] = (t[i] - a)/(b - b);
	
	for (int i = 0; i < N; i++){
		u_r[i*2+0] = u[    i  ];
		u_r[i*2+1] = u[    i+1];
		u_f[i*2+0] = u[2*N-i  ];
		u_f[i*2+1] = u[2*N-i-1];
	}
	
	// Copy solution vector	
	for (int i = 0; i < M; i++){
		sig0[i] = s[      i];
		sigM[i] = s[2*N*M+i];
	}

	for (int i = 0; i < N*M; i++){
		s_r[i] = s[        i];
		s_f[i] = s[(N+1)*M+i];
	}

	for (int i = 0; i < M; i++){
		s_m[i] = s[N*M+i];
	}
	
	// Generate boundary condition coefficients A, B, c
	for (int i = 0; i < M; i++){
		int j = 0;
		x[i] = s[i*M + j];
	}
	bcs('S', 0, N, M, a, b, p, x, A, B, c);

	/*--- CAlCULATE X FROM S (TAKE SHOTS) ---*/
	// calculate x_r (rear + midpoint)
	for (int i = 0; i < N; i++){
		double u0 = u_r[i*2+0];
		double u1 = u_r[i*2+1];
		for (int j = 0; j < M; j++){
			sig[j] = s_r[i*M+j];
		}
		// CALCULATE CHI NOW

		for (int j = 0; j < M; j++){
			x_r[i*M=j] = chi[j];
		}
		
	}

	// calculate x_f (front + midpoint)
	for (int i = 0; i < N; i++){
	}


	// START FROM HERE!!!
	// START FROM HERE!!!
	// START FROM HERE!!!
	// START FROM HERE!!!
	// START FROM HERE!!!
	// START FROM HERE!!!
	
	
	/*--- CALCULATE F (VECTOR TO BE ZEROED) ---*/
	
	
	/*--- CALCULATE DF (JACOBIAN MATRIX) ---*/



	/*--- FINALIZATION ---*/

	// free memory
	free(x    ); 
	free(dx   ); 
	
	free(u    ); 
	free(u_r  ); 
	free(u_f  ); 

	free(sig  ); 
	free(sig0 ); 
	free(sigM ); 
	free(sigp ); 
	free(s_r  ); 
	free(s_f  ); 
	free(s_m  ); 

	free(chi  ); 
	free(chip ); 
	free(x_r  ); 
	free(x_f  ); 

	free(r    ); 

	free(F    ); 
	free(F_r  ); 
	free(F_f  ); 
	free(F_m  ); 
	free(F_rc ); 

	free(DF   ); 
	free(DF_r ); 
	free(DF_f ); 
	free(DF_rc);

	free(I    ); 
	free(G    ); 
	free(A    ); 
	free(B    ); 
	free(c    ); 
}

void newt(){
}

void solve(double *s){
}

/* Given initial condition x(u0;p) = x0, integrate from u = u0 to u1
 * using a fourth-order Runge-Kutta scheme and output x(u1;p) = x1. 
 *   M : number of variables (order of the differential system)
 *   NRK : number of Runge-Kutta steps
 *   a : lower bound of domain
 *   b : upper bound of domain
 *
 * NOTE FOR LATER (5/1/17): Replace odesS0 with general odes for function evaluations.
 *
 */
void rk4(int M, int NRK, double a, double b, double u0, double u1, 
         double *p, double *x0, double *x1){
	// define pointers and allocate memory
	double *x00 = (double *) calloc(M, sizeof(double));
	double *x01 = (double *) calloc(M, sizeof(double));
	double *x02 = (double *) calloc(M, sizeof(double));
	double *f   = (double *) calloc(M, sizeof(double));

	// initialization
	double u = u0;
	for (int i = 0; i < M; i++)
		x02[i] = x0[i];
	
	// Runge-Kutta multi-step integration
	double du = (u1 - u0)/NRK;
	for (int irk = 0; irk < NRK; irk++){
		// setup
		u = irk*du;
		for (int i = 0; i < M; i++)
			x00[i] = x02[i];

		// first step
		for (int i = 0; i < M; i++)
			x01[i] = x00[i];
		odesS0(a, b, u, p, x01, f);
		for (int i = 0; i < M; i++)
			x02[i] += du*(1.0/6.0)*f[i];

		// second and third steps
		for (int j = 0; j < 2; j++){
			for (int i = 0; i < M; i++)
				x01[i] = x00[i] + du*0.5*f[i];
			odesS0(a, b, u, p, x01, f);
			for (int i = 0; i < M; i++)
				x02[i] += du*(1.0/3.0)*f[i];
		}

		// fourth step
		for (int i = 0; i < M; i++)
			x01[i] = x00[i] + du*f[i];
		odesS0(a, b, u, p, x01, f);
		for (int i = 0; i < M; i++)
			x02[i] += du*(1.0/6.0)*f[i];
	}

	// store results
	for (int i = 0; i < M; i++)
		x1[i] = x02[i];

	// free memory
	free(x00);
	free(x01);
	free(x02);
	free(f  );
}


#endif
