#include "socp/model.hpp"
#include "socp/map.hpp"

#include "opinion_params.hpp"

#ifndef _OPINION_H_
#define _OPINION_H_


/**************			opinion class			******************************/
class opinion : public model
{

public:
  opinion(){};
  virtual ~opinion(){};

  struct model_params_t model_params;
  struct problem_params_t problem_params;
  struct homotopy_params_t homotopy_params;
  struct ode_params_t ode_params;

  mstate GetStates(mstate const& X) const;
  mstate GetCostates(mstate const& X) const;
  mstate GetOpinions(mstate const& X) const;
  void Split(mstate const& X, mstate *Y, mstate *P) const;
  mstate Fuse(mstate const& Y, mstate const& P) const;

protected:
  /**
	* State model of the vehicle : dX/dt = Model(t, X)
	* @param t the time
	* @param X the state
	* @return dX/dt as a model state
	*/
	virtual mstate Model(real const& t, mstate const& X) const;

	/**
	* Control model of the vehicle : U = Control(t, X)
	* @param t the time
	* @param X the state
	* @return the control U as a control state
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	* @param t the time
	* @param X the state
	* @return the Hamiltonian value
	*/
	virtual real Hamiltonian(real const& t, mstate const& X) const;

	/**
	* Integrate state equations
	* @param X the state
	* @param t0 the initial time
	* @param tf the final time
	* @param isTrace trace flag (0 if no trace required)
	* @return the state at tf
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace);

  static mstate static_Model(real const& t, opinion::mstate const& X, void *model);
};

#endif //_OPINION_H_
