#include <math.h>
#include "mpmod.H"



mp_complex mpregkernc(mp_real x, mp_real x0, mp_real b)
{
  static mp_real pi   = mppic ;
  static mp_real pi2  = mppic*mppic;
  static mp_complex I = mp_complex(0.0,1.0);
  mp_real b2 = b*b;
  
  mp_complex y = exp(-b2*(x0-x)*(x0-x))*(cos(pi*b2*(x0-x))+I*sin(pi*b2*(x0-x)))*(I*cos(2.0*pi*b2*(x0-x))*sinh(pi2*b2)-sin(2.0*pi*b2*(x0-x))*cosh(pi2*b2));

  return y;
}

mp_complex mpcregkern(mp_real x, mp_real b)
{
  static mp_real pi   = mppic ;
  static mp_real pi2  = mppic*mppic;
  static mp_complex I = mp_complex(0.0,1.0);
  mp_real b2 = b*b;
  
  mp_complex y = exp(-b2*x*x)*(cos(pi*b2*x)+I*sin(pi*b2*x))*(I*cos(2.0*pi*b2*x)*sinh(pi2*b2)-sin(2.0*pi*b2*x)*cosh(pi2*b2));

  return y;
}

mp_real mpregkerne1(mp_real x_, mp_real b_)
{
  static mp_real pi = mppic ;
  // cout << "pi = " << pi << endl;
  mp_real _y = exp(-b_ * b_ * x_ * x_ ) * sin(pi * b_ * b_ * x_ );
  return _y;
}

mp_real mpregkerne2(mp_real x_, mp_real b_)
{
  static mp_real pi = mppic ;
  mp_real _y = exp(-b_ * b_ * x_ * x_ ) * cos(pi * b_ * b_ * x_ );
  return _y;
}

mp_real ColeCole(mp_real x_, mp_real b_, mp_real t0)
{
  static mp_real pi = mppic;
  //  mp_real factor = sin(pi*(b_))/(mp_real(2.0)*pi);
  mp_real factor = sin(pi*(b_))/(mp_real(2.0)*pi);
  mp_real x0 = -log(t0);
  mp_real y_;
  y_ = mp_real(1.0)/(cosh(b_ * (x_ - x0)) + cos(pi*b_));
  y_ *= factor;
  return y_;

}


double UpperBound(double sigma, double gamma)
{
  static double pi = 4.0 * atan(1.0);
  double bound = pow(2.0*pi,1.5) * sigma/gamma * exp(0.25/(gamma*gamma));
  // double bound = sigma/(sqrt(2.0*pi)*gamma) * exp(0.25/(gamma*gamma));
  return bound;
}
