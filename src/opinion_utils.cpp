#include "opinion_model.hpp"

opinion::mstate opinion::GetStates(mstate const& X) const
{
  mstate Y(this->model_params.N+1);

  for (int i = 0; i <= this->model_params.N; i++)
  {
    Y[i] = X[i];
  }

  return Y;
}

opinion::mstate opinion::GetCostates(mstate const& X) const
{
  mstate P(this->model_params.N+1);

  for (int i = 0; i <= this->model_params.N; i++)
  {
    P[i] = X[this->model_params.N+1+i];
  }

  return P;
}

opinion::mstate opinion::GetOpinions(mstate const& X) const
{
  mstate Y = this->GetStates(X);

  for (int i = 1; i <= this->model_params.N; i++)
  {
    Y[i] += Y[0];
  }

  return Y;
}

void opinion::Split(mstate const& X, mstate *Y, mstate *P) const
{
  *Y = this->GetStates(X);
  *P = this->GetCostates(X);
}

opinion::mstate opinion::Fuse(mstate const& Y, mstate const& P) const
{
  int N = this->model_params.N;

  mstate X(2*N+2);

  for (int i = 0; i <= N; i++)
  {
    X[i] = Y[i];
    X[i+N+1] = P[i];
  }

  return X;
}
