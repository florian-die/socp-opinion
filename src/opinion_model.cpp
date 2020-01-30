#include "opinion_model.hpp"
#include "opinion_functions.hpp"

opinion::mstate opinion::Model(real const& t, mstate const& X) const
{
  int N = this->model_params.N;

  mcontrol U = this->Control(t,X);
  real u = U[0];

  mstate Y = this->GetStates(X);

  mstate dY(N+1);
  dY[0] = u;

  for (int i = 1; i <= N; i++)
  {
    dY[i] = - opinion_functions::h(dY[i]) - u;
  }

  mstate P = this->GetCostates(X);

  mstate dP(N+1);
  dP[0] = 0.0;

  for (int i = 1; i <= N; i++)
  {
    dP[i] = P[i]*opinion_functions::dh(dY[i]);
  }

  return this->Fuse(dY,dP);
}

opinion::mcontrol opinion::Control(real const& t, mstate const& X) const
{
  mcontrol U(1);

  mstate P = this->GetCostates(X);

  U[0] = (P[1]+P[this->model_params.N]) / this->homotopy_params.u;

  return U;
}

real opinion::Hamiltonian(real const& t, mstate const& X) const
{
  mcontrol U = this->Control(t,X);
  real u = U[0];

  mstate dY = this->GetStates(this->Model(t,X));
  mstate P = this->GetCostates(X);

  real H = 1 + u*u;

  for (int i = 0; i <= this->model_params.N; i++)
  {
    H += P[i]*dY[i];
  }
}

opinion::mstate opinion::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace)
{
  real t = t0;							// time
	real dt = (tf-t0)/this->ode_params.steps;		// time step
	mstate Xs = X;

	for (int i = 0; i < this->ode_params.steps; i++)
  {
		// Solve ODE
		Xs = odeTools::RK4(t, Xs, dt, &static_Model, (void*) this);
		t += dt;
	}

	return Xs;
}

opinion::mstate opinion::static_Model(real const& t, opinion::mstate const& X, void *model)
{
	return static_cast<opinion*>(model)->Model(t, X);
}
