/*-------------------------------------------------------------*/
/*--  STEADY MOTION OF A CLOSELY FITTING VESICLE IN A TUBE:  --*/
/*--      PERTURBATIVE SOLUTION IN THE NARROW-GAP LIMIT      --*/
/*-------------------------------------------------------------*/

Description: 
  Numerical solution for the gap thickness as a function of the
	lateral coordinate by the multiple shooting method. If bending
	elasticity is neglected, then the shape equation is a third-
	order, nonlinear ODE. If bending is included, the order of the
	equation increases to five. Higher-order corrections in terms
	of the small-gap parameter are governed by linear ODEs, but
	require accurate tabulation of lower-order solutions.
	
	See supporting documentation for the derivation of the basic
	equations and solution algorithm.

Language: C++

Libraries: LAPACK, ATLAS, BLAS, GSL

Directory tree:
	include - contains header files
	src - contains input parameters and executable
	output - stores output files (in .txt format)
	postproc - contains scripts for post-processing	


EDIT EVERYTHING BELOW
EDIT EVERYTHING BELOW
EDIT EVERYTHING BELOW

Implementations (header files):
	fd.h - finite difference subroutines
	la.h - linear algebra subroutines
	rw.h - read / write subroutines

Input parameters (in src/params.in)
	CA - capillary number
	BO - bond number
	MA - Marangoni number
	TSTOP - time when sphere stops moving
	R1 - upper limit of radial domain
	T1 - upper limit of time domain 
	DR - spatial grid size
	DT - time step size
	DTREC - time points to output data file

Instructions:
	1. Make executable,
		cd src
		make clean
		make
	2. Set input parameters in src/params.in
	3. Run executable,
		cd src
		./run
	4. Post-process output files written to output.

Notes:
	Since the problem is fourth-order in space, the CFL number is given by
	
	        dt D
	  CFL = ----
	        dr^4
	
	where D is a fourth-order diffusion coefficient. Typically,
	
	  D ~ O(r1^6 / Ca)    or    D ~ O(r1^4 Ma)
	
	Having CFL = O(1) is desirable for numerical accuracy, but is not 
	required for numerical stability.
