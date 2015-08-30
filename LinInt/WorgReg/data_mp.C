#include <math.h>
#include "mpmod.H"

/* Rauschen u. Fehler */
double aint(double x)

{
    if (fabs(x) < 1.0) return 0.0;
    return floor(fabs(x))*(x >= 0 ? 1.0 : -1.0);
}

//#if	0
double randlc(double& x, const double a)

/*
!   This routine returns a uniform pseudorandom double precision number in the
!   range (0, 1) by using the linear congruential generator

!   x_{k+1} = a x_k  (mod 2^46)

!   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
!   before repeating.  The argument A is the same as 'a' in the above formula,
!   and X is the same as x_0.  A and X must be odd double precision integers
!   in the range (1, 2^46).  The returned value RANDLC is normalized to be
!   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
!   the new seed x_1, so that subsequent calls to RANDLC using the same
!   arguments will generate a continuous sequence.

!   This routine should produce the same results on any computer with at least
!   48 mantissa bits in double precision floating point data.  On Cray systems,
!   double precision should be disabled.

!   David H. Bailey     October 26, 1990
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    static int ks = 0;
    static double r23, r46, t23, t46;

/*
!   If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
!   T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
!   by merely using the ** operator, in order to insure that the results are
!   exact on all systems.  This code assumes that 0.5D0 is represented exactly.
*/
    if (ks == 0) {
	r23 = r46 = t23 = t46 = 1.0;
	for (int i = 0; i < 23; i++) {
	    r23 *= 0.5;
	    t23 *= 2.0;
	}
	for (i = 0; i < 46; i++) {
	    r46 *= 0.5;
	    t46 *= 2.0;
	}
	ks = 1;
    }

    //!   Break A into two parts such that A = 2^23 * A1 + A2 and set
    //  X = N.
    double t1 = r23*a;
    double a1 = aint(t1);
    double a2 = a-t23*a1;

    //!   Break X into two parts such that X = 2^23 * X1 + X2, compute
    //!   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    //!   X = 2^23 * Z + A2 * X2  (mod 2^46).

    t1 = r23*x;
    double x1 = aint(t1);
    double x2 = x-t23*x1;
    t1 = a1*x2+a2*x1;
    double t2 = aint(r23*t1);
    double z = t1-t23*t2;
    double t3 = t23*z+a2*x2;
    double t4 = aint(r46*t3);
    x = t3-t46*t4;
    return r46*x;
}
// #endif


double Xi1double(double x, int N = 100) // default argument: N = 100
{

  //  int N = 100;
  double xi,xk,a;
  //  double b[100],c[100];
  static double pi = 4.0 * atan(1.0);
  double T = 30.0;
  double f = 2.0 * pi/T;
  double* const b = new double [N];
  double* const c = new double [N];
  
  int j; 

  
  xk = 1.0;
  double x1 = 3.0;
  a = 3.0;
  xi = 0.0;

  for(int i = 0; i < N; i++)
    { 
      j = N  - i - 1 ;
      b[j] = randlc(xk,a);
      c[j] = randlc(x1,a);
      
    }
  
  for(i = 0; i < N; i++)
    {
      xi += b[i] * cos(i*f*x) + c[i] * sin(i*f*x);
  
      if(fabs(xi) > 10.0)
	{ xi = xi/100.0;}
      else { xi = xi;} 
    }

  if(fabs(xi) > 1.0)
    {xi = xi/10.0;}
  else 
    { xi = xi ;}
  
  delete [] b;
  delete [] c;
  return xi;
}

double Xi2double(double x, int N = 100) // default argument: N = 100
{
  //  cout << setiosflags(ios::uppercase);

  //  int N = 100;
  double xi,xk,a;
  //  double b[100],c[100];
  static double pi = 4.0 * atan(1.0);
  double T = 30.0;
  double f = 2.0 * pi/T;
  double* const b = new double [N];
  double* const c = new double [N];
  
  int j; 

  
  xk = 1.0;
  double x1 = 3.0;
  a = 3.0;
  xi = 0.0;

  for(int i = 0; i < N; i++)
    { 
      j = N  - i - 1 ;
      b[j] = randlc(xk,a);
      c[j] = randlc(x1,a);
      
    }
  
  for(i = 0; i < N; i++)
    {
      xi += b[i] * cos(i*f*x) + c[i] * sin(i*f*x);
    }

  if(fabs(xi) > 1.0)
    {xi = xi/100.0;}
  else 
    { xi = xi ;}
  
  delete [] b;
  delete [] c;
  return xi;
}

double Xi3double(double x, int N = 100) // default argument: N = 100
{
  //  cout << setiosflags(ios::uppercase);

  //  int N = 100;
  double xi,xk,a;
  //  double b[100],c[100];
  static double pi = 4.0 * atan(1.0);
  double T = 45.0;
  double f = 2.0 * pi/T;
  double* const b = new double [N];
  double* const c = new double [N];
  
  int j; 

  
  xk = 1.0;
  double x1 = 3.0;
  a = 3.0;
  xi = 0.0;

  for(int i = 0; i < N; i++)
    { 
      j = N  - i - 1 ;
      b[j] = randlc(xk,a);
      c[j] = randlc(x1,a);
      
    }
  
  for(i = 0; i < N; i++)
    {
      xi += pow(-1.0,i) * (b[i] * cos(i*f*x) + c[i] * sin(i*f*x));
    }

  if(fabs(xi) > 1.0)
    {xi = xi/10.0;}
  else 
    { xi = xi ;}
  

  delete [] b;
  delete [] c;
  return xi;
}


mp_real Xi1(mp_real x, int N = 100) // default argument: N = 100
{
  //  cout << setiosflags(ios::uppercase);

  //  int N = 100;
  //  double xi,xk,a;
  double xk,a;
  mp_real xi;
  //  double b[100],c[100];
  //  static double pi = 4.0 * atan(1.0);
  static mp_real pi = mppic;
  mp_real T = 30.0;
  mp_real f = 2.0 * pi/T;
  mp_real* const b = new mp_real [N];
  mp_real* const c = new mp_real [N];
  
  int j; 

  
  xk = 1.0;
  double x1 = 3.0;
  a = 3.0;
  xi = 0.0;

  for(int i = 0; i < N; i++)
    { 
      j = N  - i - 1 ;
      b[j] = randlc(xk,a);
      c[j] = randlc(x1,a);
      
    }
  
  for(i = 0; i < N; i++)
    {
      xi += b[i] * cos(i*f*x) + c[i] * sin(i*f*x);
  
      if(abs(xi) > 10.0)
	{ xi = xi/100.0;}
      else { xi = xi;} 
    }

  if(abs(xi) > 1.0)
    {xi = xi/10.0;}
  else 
    { xi = xi ;}
  
  delete [] b;
  delete [] c;
  return xi;
}

mp_real Xi2(mp_real x, int N = 100) // default argument: N = 100
{
  //  cout << setiosflags(ios::uppercase);

  //  int N = 100;
  double xk,a;
  mp_real xi;
  //  double b[100],c[100];
  static mp_real pi = mppic;
  mp_real T = 30.0;
  mp_real f = 2.0 * pi/T;
  mp_real* const b = new mp_real [N];
  mp_real* const c = new mp_real [N];
  
  int j; 

  
  xk = 1.0;
  double x1 = 3.0;
  a = 3.0;
  xi = 0.0;

  for(int i = 0; i < N; i++)
    { 
      j = N  - i - 1 ;
      b[j] = randlc(xk,a);
      c[j] = randlc(x1,a);
      
    }
  
  for(i = 0; i < N; i++)
    {
      xi += b[i] * cos(i*f*x) + c[i] * sin(i*f*x);
    }

  if(abs(xi) > 1.0)
    {xi = xi/100.0;}
  else 
    { xi = xi ;}
  
  delete [] b;
  delete [] c;
  return xi;
}

mp_real Xi3(mp_real x, int N = 100) // default argument: N = 100
{
  //  cout << setiosflags(ios::uppercase);

  //  int N = 100;
  double xk,a;
  mp_real xi;
  //  double b[100],c[100];
  static mp_real pi = mppic;
  mp_real T = 30.0;
  mp_real f = 2.0 * pi/T;
  mp_real* const b = new mp_real [N];
  mp_real* const c = new mp_real [N];
  
  int j; 

  
  xk = 1.0;
  double x1 = 3.0;
  a = 3.0;
  xi = 0.0;

  for(int i = 0; i < N; i++)
    { 
      j = N  - i - 1 ;
      b[j] = randlc(xk,a);
      c[j] = randlc(x1,a);
      
    }
  
  for(i = 0; i < N; i++)
    {
      xi += pow(-1.0,i)*(b[i] * cos(i*f*x) + c[i] * sin(i*f*x));
    }

  if(abs(xi) > 1.0)
    {xi = xi/10.0;}
  else 
    { xi = xi ;}
  
  delete [] b;
  delete [] c;
  return xi;
}

/* Daten, Modelle */
mp_real mpe2(mp_real x_,mp_real xi)
{ 
  // mp_real xi = mp_real(randlc(3.0,3.0));
  //  mp_real xi = Xi1(x_);;
  double x = dble(x_);
  // double xi = .5 * (Xi2double(x) + Xi3double(x));
  double nu = Xi1double(x);
  mp_real p;
  mp_real q;
  p =  exp(-x_);
  q = mp_real(1.0) + exp(- 2.0 * x_);
  mp_real y_ = p/q;
  
  
  y_ = y_ * (mp_real(1.0) + xi * nu);
  
  //y_ += mp_real(1E-02) * xi;
  
  return y_;
}

mp_real mpe2exakt(mp_real x_)
{  
  mp_real p;
  mp_real q;
  p =  exp(-x_);
  q = mp_real(1.0) + exp(- 2.0 * x_);
  mp_real y_ = p/q;
  return y_ ;
}


/* regularisierte Lösungen, ohne Datenfehler */
mp_real mpRegexakt(mp_real x_, mp_real b_)
{ 
  static mp_real SqrtPi = sqrt(mppic);
  mp_real y_ = (b_/SqrtPi) * exp( - b_ * b_ * x_ * x_);
  return y_ ;
}




