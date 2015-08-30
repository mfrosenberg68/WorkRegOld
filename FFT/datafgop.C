#include <math.h>         // Default Mathematik-Bibliothek
#include "random.H"
#include "mpmod.H"
#include "function.H"

static double  dpi = 4.0 * atan(1.0); // pi als double 


double e1data(double x)
{
  double p,q,y;
  p = 1.0;
  q = 1.0 + exp(- 2.0 * x);
  y = p/q;
  return y;
}

double e2data(double x)
{
  double p,q,y;
  p = exp(-x);
  q = 1.0 + exp(- 2.0 * x);
  y = p/q;
  return y;
}

double e1datafgop(double x, double beta, int N = 1000)
{
  double xa,xb,xi,step,nu;
  xa = -beta;
  xb =  beta;
  step = (xb - xa)/N;
  //  nu = Xi4double(x);
  nu = 0.0;
  double e1 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e1 += Fgop(xi,beta)*e1data(x-xi);
    }
  e1 *= step;
  //  e1 *= (1.0 + 0.1 * nu);

  return e1; //noch nicht normiert !!!
}

double e2datafgop(double x, double beta, int N = 1000)
{
  double xa,xb,xi,step,nu;
  xa = -beta;
  xb =  beta;
  step = (xb - xa)/N;
    //  nu = Xi4double(x);
  nu = 0.0;
  double e2 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e2 += Fgop(xi,beta)*e2data(x-xi);
    }
  e2 *= step;
  //  e2 *= (1.0 + 0.1 * nu);
  return e2; //noch nicht normiert !!!
}


double e1datafgop_exakt(double x, double beta, int N = 1000)
{
  double xa,xb,xi,step;
  xa = -beta;
  xb =  beta;
  step = (xb - xa)/N;
  double e1 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e1 += Fgop(xi,beta)*e1data(x-xi);
    }
  e1 *= step;
  return e1; //noch nicht normiert !!!
}

double e2datafgop_exakt(double x, double beta, int N = 1000)
{
  double xa,xb,xi,step;
  xa = -beta;
  xb =  beta;
  step = (xb - xa)/N;
  double e2 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e2 += Fgop(xi,beta)*e2data(x-xi);
    }
  e2 *= step;
  return e2; //noch nicht normiert !!!
}
