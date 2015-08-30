#include <math.h>         // Default Mathematik-Bibliothek
#include "random.H"
//#include "mpmod.H"
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

double e1datagauss(double x, double beta, int N = 10000)
{
  double xa,xb,xi,step,nu;
  xb = 10.0;
//  xb = 8.0/(beta*beta);
  xa =  -xb;
  step = (xb - xa)/N;
  //  nu = Xi4double(x);
  nu = 0.0;
  double e1 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e1 += GAUSS(xi,beta)*e1data(x-xi);
    }
  e1 *= step;
  //  e1 *= (1.0 + 0.01 * nu);

  return e1;
}

double e2datagauss(double x, double beta, int N = 10000)
{
  double xa,xb,xi,step,nu;
//  xb = 8.0/(beta*beta);
  xb = 10.0;
  xa =  -xb;
  step = (xb - xa)/N;
  //  nu = Xi4double(x);
  nu = 0.0;
  double e2 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e2 += GAUSS(xi,beta)*e2data(x-xi);
    }
  e2 *= step;
  //  e2 *= (1.0 + 0.01 * nu);
  return e2;
}


double e1datagauss_exakt(double x, double beta, int N = 1000)
{
  double xa,xb,xi,step;
  xb = 8.0/(beta*beta);
  xa =  -xb;
  step = (xb - xa)/N;
  double e1 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e1 += GAUSS(xi,beta)*e1data(x-xi);
    }
  e1 *= step;
  return e1; 
}

double e2datagauss_exakt(double x, double beta, int N = 1000)
{
  double xa,xb,xi,step;
  xb = 8.0/(beta*beta);
  xa =  -xb;
  step = (xb - xa)/N;
  double e2 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e2 += GAUSS(xi,beta)*e2data(x-xi);
    }
  e2 *= step;
  return e2; 
}

double e2data_multi_gauss(double x, double beta1, double t1, double beta2, double t2, double a, int N = 1000)
{
  double xa,xb,xi,step;
  //  a2 = 1.0 - a;
  //  xb = 8.0/(beta*beta);
  xb = 100.0;
  xa =  -xb;
  step = (xb - xa)/N;
  double e2 = 0.0;

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e2 += (a*GAUSS((t1 - xi),beta1) + (1.0-a) * GAUSS((t2 - xi),beta2)) *e2data(x-xi);
    }
  e2 *= step;
  return e2; 
}

double e2mod_multi_exp(double x, double beta1, double t1, double beta2, double t2, double a, int N = 1000)
{
  double xa,xb,xi,step;
  //  a2 = 1.0 - a;
  //  xb = 8.0/(beta*beta);
  xb = 100.0;
  xa =  -xb;
  step = (xb - xa)/N;
  double e2 = 0.0;
  //  t1 = -log(t1);
  //  t2 = -log(t2);

  for(int i = 0; i <= N; i++)
    {
      xi = xa + i * step;
      e2 += (a*Expmod(xi,t1,beta1) + (1.0-a) * Expmod(xi,t2,beta2)) *e2data(x-xi);
    }
  e2 *= step;
  return e2; 
}
