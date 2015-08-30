#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "mpmod.H"


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
  cout << setiosflags(ios::uppercase);

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

mp_real mpf1(mp_real x_)
{ 
  double x = dble(x_);
  double xi = Xi1double(x);
  mp_real y_ = xi;
  return y_;
}


mp_real mpf2(mp_real x_, mp_real x0_, mp_real b_)
{
 static mp_real pi = mppic ;
 mp_real _y = exp(-b_ * b_ * (x0_ -x_) *(x0_ - x_)) * cos(pi * b_ * b_ * (x0_ -x_ ));  
 return _y;
}


mp_real mpkern(mp_real x_, mp_real x0_, mp_real b_)
{
mp_real _y = mpf1(x_) * mpf2(x_, x0_, b_);
return _y;
}

mp_real mptrapez(mp_real function(mp_real _x, mp_real _x0, mp_real _a), mp_real _x0, mp_real _a, mp_real s_, mp_real h_, int n_) 
{

  mp_real sum, zu, zd;
  //  cout <<"from trapez: s = " << s_ ;
  sum = function(s_,_x0,_a);
  for (int i = 1; i <= n_ ; i++) {
    zu = s_ + i * h_ ;
    zd = s_ - i * h_ ;
    sum += function(zu,_x0,_a) + function(zd,_x0,_a);
  }
  
  return sum *= h_ ;
}

mp_real Regtheo(mp_real x_, mp_real b_)
{
  mp_real Sqrtpi = sqrt(mppic);
  mp_real y_ = b_/Sqrtpi * exp(- b_ * b_ * x_ * x_);
  return y_ ;
}

double UpperBound(double sigma, double gamma)
{
  static double pi = 4.0 * atan(1.0);
//  double bound = pow(2.0*pi,1.5) * sigma/gamma * exp(0.25/(gamma*gamma));
  double bound = .5*pow(2.0*pi,1.5) * sigma/gamma * exp(1.2337/(gamma*gamma));
//  double bound =  sigma/(pow(2.0*pi,1.5)*gamma) * exp(1.125/(gamma*gamma));
  // double bound = sigma/(sqrt(2.0*pi)*gamma) * exp(0.25/(gamma*gamma));
  return bound;
}


int main()
{
  cout << setiosflags(ios::uppercase);
  ofstream out("Xib.out");
  out  << setiosflags(ios::uppercase);

  mp_init();
  mp_real b, x0, s, h, x1, x2;
  double s_in, h_in , b_in, x0in, dfactor, dpi,x1_in, x2_in, xd;
  int n,m;
  dpi = 4.0 * atan(1.0);

  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;
  //  mp_real Sqrtpi = sqrt(mppic);
  static mp_real factor;
    
  cin >> x0in >> b_in;
  cin >> s_in >> h_in >> n;
  cin >> x1_in >> x2_in >> m;

  s = mp_real(s_in);
  h = mp_real(h_in);
  x0 = mp_real(x0in);
  b  = mp_real(b_in);
  x1 = mp_real(x1_in);
  x2 = mp_real(x2_in);

  //  factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
  //  factor *= exp(mp_real(0.25) * pi2 * b * b);
  
  dfactor = 2.0 * (b_in/pow(dpi,1.5)) * exp(b_in*b_in*dpi*dpi*.25);

  double bound, da;
  
  cout << "+++ Precision = " << mpipl << " +++" << endl;
  cout << "b  = " << b ;
  cout << "x0 = " << x0 ;
  cout << "s  = " << s;
  cout << "h  = " << h;
  cout << "x1 = " << x1;
  cout << "x2 = " << x2;
  cout << "m  = " << m << endl;
  cout << "n  = " << n << endl;
  cout << "hn = " << h * n << endl;

///  cout << "dfaktor = " << dfactor << endl;
///  cout << "Faktor  = " << factor << endl;
  //  cout << "bound   = " << bound << endl;
      
  mp_real x, h2;
    
  h2 = (x2 - x1)/mp_real(m);

  mp_real RandomReg_mp, Random_mp, Reg_mp, Regpure_mp; 
  double RandomReg, Random, Reg, Regpure; 

  for(int i = 0; i<= m; i++){

    //    x = x1 + i * h2;
    b = x1 + i * h2;
    x = x0;
    s = x0;

    xd = dble(b);
    cout << "i = " << i << '\t' << " b = " << xd << endl;

    factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
    factor *= exp(mp_real(0.25) * pi2 * b * b);
    RandomReg_mp = mptrapez(mpkern,x,b,s,h,n);
    RandomReg_mp *= factor;
    Random_mp = mpf1(x);
    mp_real diffrandom_mp = RandomReg_mp - Random_mp ;
    Regpure_mp = Regtheo(x,b);
    Reg_mp = Regpure_mp + RandomReg_mp;

    //    xd = dble(x);

    RandomReg = dble(RandomReg_mp);
    Random    = dble(Random_mp);

    Regpure = dble(Regpure_mp);
    Reg     = dble(Reg_mp);

    //cout << x << Random_mp << RandomReg_mp << diffrandom_mp << endl;
    //cout << b << Random_mp << RandomReg_mp << diffrandom_mp << endl;
    //    cout << Regpure_mp <<  Reg_mp <<  endl;
    //    cout << Regpure << " + " << RandomReg << " = " << Reg << endl;
    
    
    //(ostream&)out << xd << '\t'<< Random << '\t' << RandomReg << '\t' << Regpure << '\t' << Reg << '\t' << bound << endl; 
    
    da =1.0/(2.0*xd);
    bound   = UpperBound(.05*Random,da);

    (ostream&)out << xd << '\t'<< Random << '\t' << RandomReg  << '\t' << fabs(RandomReg) << '\t' <<  bound << endl;  
  }
  
  return 0;
}


