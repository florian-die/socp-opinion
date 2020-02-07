/*
 * opinion.hpp
 *
 *  Created on: January 30, 2018
 *      Author: Bruno HERISSE
 */

#include "socp/model.hpp"
#include "socp/map.hpp"

//#include <iostream>

#ifndef _OPINION_H_
#define _OPINION_H_

/**************			opinion class			******************************/
class opinion:public model
{

public:

	/**
	* Vehicle parameters
	*/
	struct parameters_struct{
		real u_max;				// max normalized control
		real muU;				// weight for control cost
		real muh;				// weight for interaction
	};

	/**
	* Constructor
	*/
	opinion(map *the_map, std::string the_fileTrace);

	/**
	* Destructor
	*/
	virtual ~opinion();

	/**
	* Get data pointer
	*/
	parameters_struct *GetParameterData();

	/**
	* Function value for the considered control problem
	*/
	virtual void FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const;
	virtual void FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const;

private:

	/**
	* Map
	*/
	map *my_map;

	/**
	* Vehicle data
	*/
	struct data_struct;
	data_struct *data;

	/**
	* State model of the vehicle
	*/
	virtual mstate Model(real const& t, mstate const& X) const;

	static mstate static_Model(real const& t, opinion::mstate const& X, void *model);

	/**
	* Control model of the vehicle
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	*/
	virtual real Hamiltonian(real const& t, mstate const& X) const;

	/**
	* Integrate state equations with a RK4
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace);

	/**
	* Influence of the leader
	*/
	real h0Function(real const& y) const;

	/**
	* Derivative of the influence function
	*/
	real dh0Function(real const& y) const;

	/**
	* Interactions with agents
	*/
	real hFunction(real const& y) const;

	/**
	* Derivative of the interaction function
	*/
	real dhFunction(real const& y) const;

};

#endif //_OPINION_H_
