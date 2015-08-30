#include <math.h>
//#include "mpmod.H"

static double  dpi = 4.0 * atan(1.0); // pi als double 
static double dsqrtpi = sqrt(dpi);


double GAUSS(double x, double b)
{
  double y;
  double factor = b/dsqrtpi;
  y = exp(-b*b*x*x);
  y *= factor;
  return y;
}


double UpperBound(double sigma, double gamma)
{
  static double pi = 4.0 * atan(1.0);
  double bound = pow(2.0*pi,1.5) * sigma/gamma * exp(0.25/(gamma*gamma));
  // double bound = sigma/(sqrt(2.0*pi)*gamma) * exp(0.25/(gamma*gamma));
  return bound;
}
