#ifndef __OPINION_PARAMS__
#define __OPINION_PARAMS__

struct model_params_t
{
  int N; // number of agents
  double sigma; // max norm of control
};

struct problem_params_t
{
  double eta; // final neighborhood
};

struct homotopy_params_t
{
  double u;
  double h;
};

struct ode_params_t
{
  int steps;
};

#endif
