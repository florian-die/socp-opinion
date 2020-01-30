#ifndef __OPINION_FUN__
#define __OPINION_FUN__

namespace opinion_functions
{
  double f(double x)
  {
    return 1.0/(1.0+x*x);
  }

  double df(double x)
  {
    return -2.0*x*f(x)*f(x);
  }

  double ddf(double x)
  {
    double fx = f(x);
    return 8.0*x*x*fx*fx*fx-2.0*fx*fx;
  }

  double h(double x)
  {
    return x*f(x);
  }

  double dh(double x)
  {
    return f(x)+x*df(x);
  }

  double ddh(double x)
  {
    return 2*df(x)+x*ddf(x);
  }
}

#endif
