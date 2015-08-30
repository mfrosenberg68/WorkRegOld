#include <math.h>         // Default Mathematik-Bibliothek
#include "complex.H"      
#include "random.H"

const DComplex I = DComplex(0.0,1.0);
static double  dpi = 4.0 * atan(1.0); // pi als double 

DComplex edata(double x, double beta, double tau0)
{
  double xi;
  //  xi = Xi1double(x);
  xi = 0.0;
  DComplex g;
  //  g = 1.0/(1.0 - power(I * tau0 * exp(-x),beta));
  g = 1.0/(1.0 + exp(I*0.5*dpi*beta) * exp(-beta*x)* power(DComplex(tau0,0.0),beta));
  //  g = g + 1E-02 * xi *(real(g) +I *imag(g));
  return g;

}

DComplex edata_exakt(double x, double beta, double tau0)
{
  DComplex g;
  //  g = 1.0/(1.0 - pow(tau0 * exp(-x),beta)*(cos(0.5*dpi*beta)+I*sin(0.5*dpi*beta )));
  g = 1.0/(1.0 + exp(I*0.5*dpi*beta) * exp(-beta*x)* power(DComplex(tau0,0.0),beta));
  //  g = 1.0/(1.0 - power(I * tau0 * exp(-x),beta));
  return g;
}

double ColeCole(double x_, double b_, double t0)
{
  static double pi =  4.0 * atan(1.0);
  //  double factor = sin(pi*(b_))/(double(2.0)*pi);
  double factor = sin(pi*(b_))/((2.0)*pi);
  double x0 = -log(t0);
  double y_;
  y_ = 1.0/(cosh(b_ * (x_ - x0)) + cos(pi*b_));
  y_ *= factor;
  return y_;

}

