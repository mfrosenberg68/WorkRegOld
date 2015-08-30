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
  g = 1.0/power(1.0 - I * tau0 * exp(-x),beta);
  g = g + 1E-02 * xi *(real(g) +I *imag(g));
  return g;

}

DComplex edata_exakt(double x, double beta, double tau0)
{
  DComplex g;
  g = 1.0/power(1.0 - I * tau0 * exp(-x),beta);
  return g;
}
