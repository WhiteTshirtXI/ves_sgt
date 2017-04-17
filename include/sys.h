/* SYSTEM
 *  Class that specifies the boundary-value problem
 */

#ifndef SYS_H
#define SYS_H

class Sys {
	public:
		Sys();
		~Sys();
	
	char id[2];		// identifier that specifies model assumptions
	int order;		// order of the perturbation equations
};

// Constructor
Sys::Sys(){
}

// Destructor
Sys::~Sys(){
}

#endif
