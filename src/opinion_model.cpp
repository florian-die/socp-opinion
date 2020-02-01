#include "opinion_model.hpp"
#include "opinion_functions.hpp"

#include <cmath>

opinion::mstate opinion::Model(real const& t, mstate const& X) const
{
  int N = this->model_params.N;

  // control computation
  mcontrol U = this->Control(t,X);
  real u = U[0];

  // opinions dynamics
  mstate Y = this->GetStates(X);

  mstate dY(N+1);
  dY[0] = u;

  for (int i = 1; i <= N; i++)
  {
    dY[i] = - opinion_functions::h(Y[i]) - u;
  }

  // costates dynamics
  mstate P = this->GetCostates(X);

  mstate dP(N+1);
  dP[0] = 0.0;

  for (int i = 1; i <= N; i++)
  {
    dP[i] = P[i]*opinion_functions::dh(Y[i]);
  }

  return this->Fuse(dY,dP);
}

opinion::mcontrol opinion::Control(real const& t, mstate const& X) const
{
  mcontrol U(1);

  mstate P = this->GetCostates(X);

  real u = (P[1]+P[this->model_params.N]) / this->homotopy_params.u;

  double sigma = this->model_params.sigma;

  // saturation
  if (u > sigma || u < -sigma)
  {
    u = sigma * u / fabs(u);
  }

  U[0] = u;

  return U;
}

real opinion::Hamiltonian(real const& t, mstate const& X) const
{
  mcontrol U = this->Control(t,X);
  real u = U[0];

  mstate dY = this->GetStates(this->Model(t,X));
  mstate P = this->GetCostates(X);

  int N = this->model_params.N;

  // TODO : cost function ?
  real H = 1.0 + u*u/this->homotopy_params.u/2.0;

  for (int i = 0; i <= this->model_params.N; i++)
  {
    H += P[i]*dY[i];
  }

  // mstate Y = this->GetStates(X);
  // H += - u*(P[1]+P[N]);
  // H += - P[1]*opinion_functions::h(Y[1]);
  // H += - P[N]*opinion_functions::h(Y[N]);

  return H;
}

opinion::mstate opinion::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace)
{
  real t = t0; // time
	real dt = (tf-t0)/this->ode_params.steps;		// time step
	mstate Xs = X;

  if (isTrace)
  {
    this->output_params.file_stream.open(this->output_params.file_name.c_str(), std::ios::trunc);
    Trace(t, Xs, this->output_params.file_stream);
  }

	for (int i = 0; i < this->ode_params.steps; i++)
  {
		// Solve ODE
		Xs = odeTools::RK4(t, Xs, dt, &static_Model, (void*) this);
		t += dt;

    if (isTrace)
    {
      Trace(t, Xs, this->output_params.file_stream);
    }
	}

  if (isTrace)
  {
    this->output_params.file_stream.close();
  }

	return Xs;
}

opinion::mstate opinion::static_Model(real const& t, opinion::mstate const& X, void *model)
{
	return static_cast<opinion*>(model)->Model(t, X);
}
