/* FUNCTION EVALUATION
 *  Evaluate right-hand side vector f in the system of ODEs dy/dt = f(y,t).
 *
 * REFERENCES
 *  Secomb et al, J Fluid Mech (1986)
 *  
 * PARAMETERS
 *  ny					[input]		number of ODEs
 *  t						[input]		abscissas
 *  y						[input]		solution
 *  f						[input]		function
 *  r   = y[0]  [input]		radius
 *  psi = y[1]  [input]		tilt angle
 *  cs  = y[2]  [input]		meridional curvature
 *  qs  = y[3]  [input]		transverse shear tension
 *  p   = y[4]  [input]		pressure
 *  sig = y[5]  [input]		mean tension
 *  A   = y[6]  [input]		total surface area
 *  V   = y[7]  [input]		total volume
 *  Q   = y[8]  [input]		leakback flow rate
 *  S   = y[9]  [input]		total meridional arc length
 */

#ifndef FUNC_H
#define FUNC_H


/* HEADER FILES */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */
void func(int ny, double Ca,
          double t, double *y, double *f){
  int i;
  
  if (ny != 10)
    cout << "Error: ny should equal 10." << endl;
  
  // define variables
  double r   = y[0];
  double psi = y[1];
  double cs  = y[2];
  double qs  = y[3];
  double p   = y[4];
  double sig = y[5];
  double A   = y[6];
  double V   = y[7];
  double Q   = y[8];
  double S   = y[9];
	
	// check bounds on radius
	if (r > 0.99999999)
		r = 0.99999999;
	
	if (r <= 0)
		r = 1e-12;
	
	/*-------- RESULTS FROM LUBRICATION THEORY --------*/
	
	double Dr   = 1.0 - r  ; 
  double Dr2  = Dr  *Dr  ; 
  double Dr3  = Dr  *Dr2 ;
  double Dr4  = Dr2 *Dr2 ;
  double Dr5  = Dr2 *Dr3 ;
  double Dr6  = Dr3 *Dr3 ;
  double Dr7  = Dr3 *Dr4 ;
  double Dr8  = Dr4 *Dr4 ;
  double Dr9  = Dr4 *Dr5 ;
  double Dr10 = Dr5 *Dr5 ;
  double Dr11 = Dr5 *Dr6 ;
  double Dr12 = Dr6 *Dr6 ;
  double Dr13 = Dr6 *Dr7 ;
  double Dr14 = Dr7 *Dr7 ;
  double Dr15 = Dr7 *Dr8 ;
  double Dr16 = Dr8 *Dr8 ;
  double Dr17 = Dr8 *Dr9 ;
  double Dr18 = Dr9 *Dr9 ;
  double Dr19 = Dr9 *Dr10;
  double Dr20 = Dr10*Dr10;
  	
	double r2  = r*r;
  double Qpi = Q/M_PI;
  double g, e;
  
	// calculate local axial pressure gradient, g = (dp/dx)/Ca
	if (r > 0.02){
  	double log = gsl_sf_log(r  );

		g  = -(16.0/(1.0 - r2));
		g *= (0.5 + (1.0 - r2)/(4.0*log) + 0.5*Qpi);
		g /= (1.0 + r2 + (1.0 - r2)/log);
	}
	else {
		g  = -8.0*(1.0 + Qpi); // limiting value as r --> 0
	}

	// calculate local membrane shear rate, e = du/dr
//	if (r > 0.02){
//  	double log = gsl_sf_log(r);
//		double rlog = r*log;
//		double irlog = 1.0/rlog;
//
//		e = 0.25*g*((1 - r2)*irlog + 2.0*r) + irlog;
//	}
//	else {
  	e  = (3.0/Dr2)*Qpi                 ;
  	e += (2.0/Dr )*(1.0 + Qpi)         ;
  	e +=       1.0      + 1.45    *Qpi ;
  	e += Dr  *(0.7      + 1.15    *Qpi);
  	e += Dr2 *(0.55     + 0.981786*Qpi);
  	e += Dr3 *(0.469286 + 0.882857*Qpi);
  	e += Dr4 *(0.423214 + 0.820839*Qpi);
//		e += Dr5 *(0.395036 + 0.778982*Qpi); // higher order terms just
//		e += Dr6 *(0.376375 + 0.748554*Qpi); // emphasize the singular
//		e += Dr7 *(0.362969 + 0.724934*Qpi); // nature of the solution
//		e += Dr8 *(0.352610 + 0.705622*Qpi); // near the poles, which is
//		e += Dr9 *(0.344131 + 0.689218*Qpi); // unphysical ...
//		e += Dr10*(0.336895 + 0.674906*Qpi);
//		e += Dr11*(0.330541 + 0.662184*Qpi);
//		e += Dr12*(0.324853 + 0.650726*Qpi);
//		e += Dr13*(0.319695 + 0.640305*Qpi);
//		e += Dr14*(0.314973 + 0.630758*Qpi);
//		e += Dr15*(0.310620 + 0.621958*Qpi);
//		e += Dr16*(0.306587 + 0.613808*Qpi);
//		e += Dr17*(0.302833 + 0.606228*Qpi);
//		e += Dr18*(0.299326 + 0.599152*Qpi);
//		e += Dr19*(0.296039 + 0.592525*Qpi);
//		e += Dr20*(0.292949 + 0.586299*Qpi);
//	}
	
	/*-------------------------------------------------*/
	
	if (r > 0.02){ // check if far from end caps
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
  	double cosr = cos/r;
  	double sinr = sin/r;

  	// calculate function
  	f[0] = -sin;
  	f[1] = -cs;
  	f[2] = sinr*cosr
  	     + cs*sinr 
  	     - qs;
  	f[3] = p
  	     - sig*(cs - cosr)
  	     + 0.5*(cs + cosr)*(cs*cs - cosr*cosr)
  	     + qs*sinr;
  	f[4] = Ca*g*cos;
  	f[5] = -Ca*e;
  	f[6] = 2*M_PI*r;
  	f[7] = M_PI*r2*cos;
		f[8] = 0;
  	f[9] = 0;
	}
	else { // approximate cs, cphi as constant (spherical cap)
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
		
		f[0] = -sin;
		f[1] = -cs;
		f[2] = 0;
		f[3] = 0.5*(p - 2.0*cs*sig);
		f[4] = Ca*g*cos; 
		f[5] = -Ca*e;
		f[6] = 2.0*M_PI*r;
		f[7] = M_PI*r*r*cos;
		f[8] = 0; 
		f[9] = 0;
	}

	for (i = 0; i < ny; i++){
		f[i] *= S;
	}
}



#endif
