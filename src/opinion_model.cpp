#include "opinion_model.hpp"
#include "opinion_functions.hpp"

#include <cmath>

opinion::mstate opinion::Model(real const& t, mstate const& X) const
{
  using namespace opinion_functions;

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
    // without interaction term
    dY[i] = - h(Y[i]) - u;

    // interactions terms
    if (this->homotopy_params.h > 0.0)
    {
      for(int k = 0; k <= N; k++)
      {
        dY[i] += h(Y[k]-Y[i]) * this->homotopy_params.h;
      }
    }
  }

  // costates dynamics
  mstate P = this->GetCostates(X);

  mstate dP(N+1);
  dP[0] = 0.0;

  for (int i = 1; i <= N; i++)
  {
    // without interaction term
    dP[i] = P[i]*dh(Y[i]);

    // interactions terms
    if (this->homotopy_params.h > 0.0)
    {
      for(int k = 0; k <= N; k++)
      {
        dP[i] += P[i]*dh(Y[k]-Y[i]) * this->homotopy_params.h;

        if(i != k)
        dP[i] += - P[k]*dh(Y[k]-Y[i]) * this->homotopy_params.h;
      }
    }
  }

  return this->Fuse(dY,dP);
}

opinion::mcontrol opinion::Control(real const& t, mstate const& X) const
{
  mcontrol U(1);

  mstate P = this->GetCostates(X);

  real u = 0.0;

  // control for quadric control norm cost function
  if (this->homotopy_params.u > 0.0)
  {
    u = (P[1]+P[this->model_params.N]) / this->homotopy_params.u;
  }

  // control for time optimal formulation without interaction
  if (this->homotopy_params.u == 0.0 && this->homotopy_params.h == 0.0)
  {
    double t1 = this->solution_params.switching_times[0];

    if (t <= t1)
    {
      u = this->model_params.sigma;
    }

    if (t > t1)
    {
      u = this->SingularControl(X);
    }
  }

  // saturation
  double sigma = this->model_params.sigma;
  if (u > sigma || u < -sigma)
  {
    u = sigma * u / fabs(u);
  }

  U[0] = u;

  return U;
}

real opinion::SingularControl(mstate const& X) const
{
  int N = this->model_params.N;

  mstate Y = this->GetStates(X);

  using namespace opinion_functions;

  real u = ddh(Y[1])*h(Y[1]) - ddh(Y[N])*h(Y[N]);
  u = u / (ddh(Y[N]) - ddh(Y[1]));

  return u;
}

real opinion::Hamiltonian(real const& t, mstate const& X) const
{
  mcontrol U = this->Control(t,X);
  real u = U[0];

  mstate dY = this->GetStates(this->Model(t,X));
  mstate P = this->GetCostates(X);

  int N = this->model_params.N;

  // integral cost function part
  real H = 1.0 + u*u*this->homotopy_params.u/2.0;

  // system dynamics part
  for (int i = 0; i <= this->model_params.N; i++)
  {
    H += P[i]*dY[i];
  }

  return H;
}

opinion::mstate opinion::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace)
{
  real t = t0; // time
	real dt = (tf-t0)/this->ode_params.steps;		// time step
	mstate Xs = X;

  if (isTrace)
  {
    this->output_params.file_stream.open(this->output_params.file_name.c_str(), std::ios::app);
    this->Trace(t, Xs, this->output_params.file_stream);
    // this->output_params.file_stream.close();
  }

	for (int i = 0; i < this->ode_params.steps; i++)
  {
		// Solve ODE
		Xs = odeTools::RK4(t, Xs, dt, &static_Model, (void*) this);
		t += dt;

    if (isTrace)
    {
      // this->output_params.file_stream.open(this->output_params.file_name.c_str(), std::ios::app);
      this->Trace(t, Xs, this->output_params.file_stream);
      // this->output_params.file_stream.close();
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

void opinion::SwitchingTimesUpdate(std::vector<real> const& switchingTimes)
{
  // this->solution_params.switching_times.resize(switchingTimes.size());
	// for (int i = 0; i<switchingTimes.size(); i++)
  // {
  //   this->solution_params.switching_times[i] = switchingTimes[i];
  // }


  if (!switchingTimes.empty())
  {
    this->solution_params.switching_times[0] = switchingTimes[0];
  }
}

void opinion::SwitchingTimesFunction(real const& t, mstate const& X, real& fvec) const
{
  mstate P = this->GetCostates(X);

  fvec = P[1] + P[this->model_params.N];
}
