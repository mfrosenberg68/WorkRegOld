#include <math.h>
//#include "mpmod.H"

static double  dpi = 4.0 * atan(1.0); // pi als double 
static double dsqrtpi = sqrt(dpi);

/* Friedrichscher Operator, NICHT NORMIERT!!! */
double Fgop(double x, double b, double epsilon = 1E-10) 
{
  double t, y, t2, xb_;
  xb_ = x/b;
  t = 1.0 - (xb_ * xb_) ;
  /* 
  if(fabs(t) <= epsilon)
  {
    t2 = t;
    t = 0.0;
    cout << "+++ from fgop: xb = "<< xb_ << endl;
    cout << "+++ from fgop: t2 = "<< t2 << endl;
    cout << "+++ from fgop: t  = "<< t << endl ; 
    cout << '\v';
  }
  */
  if(t > 0.0) 
    {
      y = exp(- 1.0/t) ;
    }
  if(t <= 0.0) y = 0.0;
  
  return y;

}

double Norm_fgop(int n = 100000)
{
  double c, zu, zd;
  double h = 2.0/n;
  //  cout << "from Norm_fgop; n = " << n << endl;

  c = Fgop(0.0,1.0);
  for(int i = 1; i <= n; i++)
    {
      zu =   i * h;
      zd = - i * h;
      c  += Fgop(zd,1.0) + Fgop(zu,1.0);
    }
  c *=h;
  //  cout << "from Norm_fgop: c0 = " << c << endl;
  //  cout << '\v';
  
  return c;
}


double GAUSS(double x, double b)
{
  double y;
  double factor = b/dsqrtpi;
  y = exp(-b*b*x*x);
  y *= factor;
  return y;
}

double Expmod(double x, double t, double b)
{
  //  double arg= exp(-(fabs(x-t)));
  return b*exp(- b*fabs(exp(x)-t))*exp(x);
}


double UpperBound(double sigma, double gamma)
{
  static double pi = 4.0 * atan(1.0);
  double bound = pow(2.0*pi,1.5) * sigma/gamma * exp(0.25/(gamma*gamma));
  // double bound = sigma/(sqrt(2.0*pi)*gamma) * exp(0.25/(gamma*gamma));
  return bound;
}
