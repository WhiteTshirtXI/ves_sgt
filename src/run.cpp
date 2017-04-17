/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <math.h>
#include <string>
#include "../include/odes.h"
#include "../include/bcs.h"
#include "../include/sys.h"

using namespace std;

int main(){
	int M = 14;

	// allocate memory
	double *p   = (double*) calloc(1, sizeof(double));
	double *x   = (double*) calloc(M, sizeof(double));
	double *dx  = (double*) calloc(M, sizeof(double));

	// domain
	double t, t0, t1;
	t0 = 0.0;
	t1 = 1.0;

	// specify the BVP
	Sys bvp;	
	bvp.id[0] = 'S';
	bvp.order = 0;


	odesS0(t0, t1, t, p, x, dx);
	
	// free memory
	free(p );
	free(x );
	free(dx);
	
	return(0);
}
