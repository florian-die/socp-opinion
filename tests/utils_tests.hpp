#include "UnitTest++/UnitTest++.h"

#include "opinion_model.hpp"

SUITE(OpinionUtils)
{
  class UtilsFixture
  {
  public:

    opinion my_opinion;
    const int N = 2;
    opinion::mstate X;

    UtilsFixture() : my_opinion(this->N)
    {
      this->my_opinion.model_params.N = this->N;

      this->X = opinion::mstate{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    }
  };

  TEST_FIXTURE(UtilsFixture,GetStates)
  {
    opinion::mstate Y = this->my_opinion.GetStates(this->X);

    for (int i = 0; i <= this->N; i++)
    {
      CHECK_EQUAL(1.0,Y[i]);
    }
  }

  TEST_FIXTURE(UtilsFixture,GetCostates)
  {
    opinion::mstate P = this->my_opinion.GetCostates(this->X);

    for (int i = 0; i <= this->N; i++)
    {
      CHECK_EQUAL(2.0,P[i]);
    }
  }

  TEST_FIXTURE(UtilsFixture,GetOpinions)
  {
    opinion::mstate Y = this->my_opinion.GetOpinions(this->X);

    for (int i = 1; i <= this->N; i++)
    {
      CHECK_EQUAL(2.0,Y[i]);
    }
  }

  TEST_FIXTURE(UtilsFixture,Fuse)
  {
    opinion::mstate Y = this->my_opinion.GetStates(this->X);
    opinion::mstate P = this->my_opinion.GetCostates(this->X);

    opinion::mstate Z = this->my_opinion.Fuse(Y,P);

    for (int i = 0; i <= this->X.size(); i++)
    {
      CHECK_EQUAL(X[i],Z[i]);
    }
  }
}
