/*
 * opinion.cpp
 *
 *  Created on: January 30, 2018
 *      Author: Bruno HERISSE
 */
#include <fstream>

#include <math.h>

#include "opinion.hpp"

/**
* Vehicle data
*/
struct opinion::data_struct{
	int n;							// state dimension
	real pi;						// pi
	opinion::parameters_struct *parameters;	// those parameters can be used for continuation
	int  stepNbr;					// step number for ModelInt
	std::string strFileTrace;				// trace file
	std::ofstream fileTrace;				// trace file stream
};

/**
* Constructor
*/
opinion::opinion(map *the_map, std::string the_fileTrace) : model(7) {

	/// new vehicle
	my_map = the_map;					// map

	/// vehicle data
	data = new data_struct;
	data->n = dim;						// state dimension
	data->pi = M_PI;					// pi
	data->parameters = new parameters_struct;
		data->parameters->u_max = 5;				// max normalized control
		data->parameters->muU = 1;				// weight for time cost
		data->parameters->muh = 0.1;				// weight for interaction
	data->stepNbr = 100;				// step number for ModelInt default = 2000
	data->strFileTrace = the_fileTrace;	// trace file

	// trace file
  	data->fileTrace.open(data->strFileTrace.c_str(), std::ios::trunc);	// erase file
};

/**
* Destructor
*/
opinion::~opinion(){

	// close trace file
  	data->fileTrace.close();

	delete(data);
	delete(my_map);
};

/**
* State model of the vehicle
*/
opinion::mstate opinion::Model(real const& t, mstate const& X) const{
	mstate Xdot(X.size());

    // control computation
	real sumP = 0;
	for (int k=1;k<data->n; k++){
		sumP += X[k+data->n];
	}
	real u = sumP / data->parameters->muU;
	real norm_u = sqrt(u*u);
	if (norm_u > data->parameters->u_max){
		u = u/norm_u*data->parameters->u_max;
		norm_u = data->parameters->u_max;
	}

	// state and costate equations
    Xdot[0] = u;
	for (int k=1;k<data->n; k++){
		real sumh = 0;
		for (int j=1;j<data->n; j++){
			if (j!=k) sumh += hFunction(X[j]-X[k]);
		}
		Xdot[k] = sumh - h0Function(X[k]) - u;
	}
	Xdot[data->n] = 0;
	for (int k=1;k<data->n; k++){
		real sumdh = 0.0;
		for (int j=1;j<data->n; j++){
			if (j!=k) sumdh += (X[k+data->n] - X[j+data->n])*dhFunction(X[j]-X[k]); // PAS D'ERREUR !
		}
		Xdot[k+data->n] = sumdh + X[k+data->n]*dh0Function(X[k]);
	}

	return Xdot;

};

opinion::mstate opinion::static_Model(real const& t, opinion::mstate const& X, void *model){ // static model to be passed as an argument of the ODE solver
	return static_cast<opinion*>(model)->Model(t, X);
}

/**
* Control model of the vehicle
*/
opinion::mcontrol opinion::Control(real const& t, mstate const& X) const{
    // control computation
	real sumP = 0;
	for (int k=1;k<data->n; k++){
		sumP += X[k+data->n];
	}
	real u = sumP / data->parameters->muU;
	real norm_u = sqrt(u*u);
	if (norm_u>data->parameters->u_max){
		u = u/norm_u*data->parameters->u_max;
		norm_u = data->parameters->u_max;
	}

	mcontrol control(1);
	control[0] = u;

	return control;

}

/**
* Hamiltonian of the vehicle
*/
real opinion::Hamiltonian(real const& t, mstate const& X) const{
	// control computation
	real sumP = 0;
	for (int k=1;k<data->n; k++){
		sumP += X[k+data->n];
	}
	real u = sumP / data->parameters->muU;
	real norm_u = sqrt(u*u);
	if (norm_u>data->parameters->u_max){
		u = u/norm_u*data->parameters->u_max;
		norm_u = data->parameters->u_max;
	}

	// H is computed
	real H = 1 + data->parameters->muU*norm_u*norm_u/2 - sumP*u;
	for (int k=1;k<data->n; k++){
		real sumh = 0;
		for (int j=1;j<data->n; j++){
			if (j!=k) sumh += hFunction(X[j]-X[k]);
		}
		H += X[k+data->n]*(sumh-h0Function(X[k]));
	}

	return H;
}

/**
* Integrate state equations with an ODE solver
*/
opinion::mstate opinion::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace){
	real t = t0;							// time
	real dt = (tf-t0)/data->stepNbr;		// time step
	mstate Xs = X;

	if (isTrace) Trace(t, Xs, data->fileTrace);

  	for (int i = 0; i < data->stepNbr; i++) {
		// Solve ODE
		Xs = odeTools::RK4(t, Xs, dt, &static_Model, (void*) this);
		t += dt;

		if (isTrace) Trace(t, Xs, data->fileTrace);
	}

	return Xs;
};

/**
* Get parameters pointer
*/
opinion::parameters_struct *opinion::GetParameterData(){
	return data->parameters;
}

/**
* Function value for the considered control problem
*/
void opinion::FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const{
	// Compute the function to solve
	for (int j=0;j<data->n;j++){
		if (mode_X[j]==1){
			// continuity => X(k+n) = 0 (transversality condition)
			fvec[j] = X_tf[j+data->n];
		}else{
			// continuity => X(k) = Xf(k)
			fvec[j] = X_tf[j] - Xf[j];
		}
	}
};

void opinion::FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const{
	// Compute the function to solve
	for (int j=0;j<data->n;j++){
		if (mode_X[j]==1){
			// continuity => X(k+n) = 0 (transversality condition)
			fvec[j] = X_tf[j+data->n];
		}else{
			// continuity => X(k) = Xf(k)
			fvec[j] = X_tf[j] - Xf[j];
		}
	}
	fvec[data->n] = Hamiltonian(tf,X_tf);
};

/**
* Influence of the leader
*/
real opinion::h0Function(real const& y) const{
	return y/(1+y*y);
}

/**
* Derivative of the influence function
*/
real opinion::dh0Function(real const& y) const{
	return (1-y*y)/(1+y*y)/(1+y*y);
}

/**
* Interactions with agents
*/
real opinion::hFunction(real const& y) const{
	return data->parameters->muh*y/(1+y*y);
}

/**
* Derivative of the interaction function
*/
real opinion::dhFunction(real const& y) const{
	return data->parameters->muh*(1-y*y)/(1+y*y)/(1+y*y);
}
