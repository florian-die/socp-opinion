#ifndef __OPINION_PARAMS__
#define __OPINION_PARAMS__

#include <vector>

struct model_params_t
{
  int N; // number of agents
  double sigma; // max norm of control
};

struct problem_params_t
{
  double eta; // final neighborhood
};

struct solution_params_t
{
    std::vector<double> switching_times;
};

struct homotopy_params_t
{
  double u;
  double h;
};

struct ode_params_t
{
  int steps;
  std::string output_file;
};

struct output_params_t
{
  std::string file_name;
  std::ofstream file_stream;
};

struct params_t
{
  struct model_params_t model;
  struct problem_params_t problem;
  struct solution_params_t solution;
  struct homotopy_params_t homotopy;
  struct ode_params_t ode;
  struct output_params_t output;
};

#endif
