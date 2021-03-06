#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "mpmod.H"
#include "tarray.h"       //DOLIB DoubleArrayx :-)
#include "vtype.h"
#include "qc_utilities.H"


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

double Xi2(double x, int N = 100) // default argument: N = 100
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
    }

  if(fabs(xi) > 1.0)
    {xi = xi/100.0;}
  else 
    { xi = xi ;}
  
  delete [] b;
  delete [] c;
  return xi;
}

double Xi3(double x, int N = 100) // default argument: N = 100
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
    }


  if(fabs(xi) > 1.0)
    {xi = xi/10.0;}
  else 
    { xi = xi ;}

  

  delete [] b;
  delete [] c;
  return xi;
}


double Xi4(double x, int N = 100) // default argument: N = 100
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


mp_real mpe2(mp_real x_)
{
  mp_real p;
  mp_real q;
  p =  exp(- x_) ;
  q = mp_real(1.0) + exp(- 2.0 * x_);
  mp_real y_ = p/q;
  return y_;
}

double e2(double x_)
{
  double y_ = exp(-x_)/(1.0 + exp(-2.0 * x_));
  return y_ ;
}

double kerngauss(double x_ , double x0_, double g)
{
  double y_ = e2(x_) * exp(- g * g * (x_ - x0_) * (x_ - x0_));
  return y_ ;
}

double trapez(double function(double _x, double _x0, double _a), double _x0, double _a, double s_, double h_, int n_) 
{

  double sum, zu, zd;
  //  cout <<"from trapez: s = " << s_ ;
  sum = function(s_,_x0,_a);
  for (int i = 1; i <= n_ ; i++) {
    zu = s_ + i * h_ ;
    zd = s_ - i * h_ ;
    sum += function(zu,_x0,_a) + function(zd,_x0,_a);
  }
  
  return sum *= h_ ;
}



mp_real mpe2gauss(mp_real x0_, mp_real g)
{ 
  static mp_real SQRTpi = sqrt(mppic) ;
  //double gd = 1.0;  /* gd = 1/2 gamma !!!; hier gamma = 0.5 */
  double x0d = dble(x0_);
//  double xi = Xi4(x0d); /*Rauschen*/
  double xi = 0.0; /*Rauschen*/
  double gd  = dble(g);
  mp_real fac = g/SQRTpi;
  double h = 1E-04;
  int n = 100000;
  double yd_ = trapez(kerngauss,x0d,gd,x0d,h,n); 
  mp_real y_ = fac * mp_real(yd_);
  //y_ = y_ * (mp_real(1.0) + mp_real(0.0) * xi);
  return y_;
}

mp_real mpe2gausstheo(mp_real x0_, mp_real g)
{ 
  static mp_real SQRTpi = sqrt(mppic) ;
//  double gd = 1.0;  /* gd = 1/2 gamma !!!; hier gamma = 0.5 */
  double x0d = dble(x0_);
  double gd  = dble(g);
  mp_real fac = g/SQRTpi;
  double h = 1E-04;
  int n = 100000;
  double yd_ = trapez(kerngauss,x0d,gd,x0d,h,n); 
  mp_real y_ = fac * mp_real(yd_);
  return y_;
}


mp_real mpf2(mp_real x_, mp_real x0_, mp_real b_)
{
 static mp_real pi = mppic ;
 mp_real _y = exp(-b_ * b_ * (x0_ -x_) *(x0_ - x_)) * cos(pi * b_ * b_ * (x0_ -x_ ));  
 return _y;
}


mp_real mpkern(mp_real x_, mp_real x0_, mp_real b_, mp_real g_)
{
mp_real _y = mpe2gauss(x_,g_) * mpf2(x_, x0_, b_);
return _y;
}

mp_real mptrapez(mp_real function(mp_real _x, mp_real _x0, mp_real _a, mp_real g), mp_real _x0, mp_real _a, mp_real g, mp_real s_, mp_real h_, int n_) 
{

  mp_real sum, zu, zd;
  //  cout <<"from trapez: s = " << s_ ;
  sum = function(s_,_x0,_a,g);
  for (int i = 1; i <= n_ ; i++) {
    zu = s_ + i * h_ ;
    zd = s_ - i * h_ ;
    sum += function(zu,_x0,_a,g) + function(zd,_x0,_a,g);
  }
  
  return sum *= h_ ;
}


mp_real Regtheo(mp_real x_, mp_real b_, mp_real g)
{
  mp_real Sqrtpi = sqrt(mppic);
  mp_real a2 = (b_ * b_) + (g * g) ; // a^2 + g^2
  mp_real a  = sqrt(a2);
  mp_real gb =  (g*b_)/a;
  mp_real y_ = gb/Sqrtpi * exp(- gb * gb * x_ * x_);
  return y_ ;
}

mp_real GaussDichte(mp_real x_,mp_real g)
{
  mp_real Sqrtpi = sqrt(mppic);
  mp_real y_ = exp(- g * g * x_ * x_);
  y_ *= g/Sqrtpi;
  return y_ ;
}

int main()
{
   
  cout  << setiosflags(ios::uppercase);
  /*
  ofstream out("GAUSS.out");
  out   << setiosflags(ios::uppercase);
  ofstream dout("GAUSSMOD.out");
  dout  << setiosflags(ios::uppercase);
  */

  mp_init();
  double s_in, h_in , b_in, g_in, dfactor, dpi, x1_in, x2_in;
  int N, m, Control;

  dpi = 4.0 * atan(1.0);
  
  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;
  static mp_real factor;
  
  mp_real a, b, g, x0, s, h, x1, x2;

  /* Parameter als doubles einlesen */
  cin >> b_in >> g_in >> Control;
  cin >> s_in >> h_in >> N;
  cin >> x1_in >> x2_in >> m;

  /* Konversion: double -> mp_real */
  b  = mp_real(b_in);
  g  = mp_real(g_in);
  s  = mp_real(s_in);
  h  = mp_real(h_in);
  x1 = mp_real(x1_in);
  x2 = mp_real(x2_in);
  a  = 1/(mp_real(2.0)*b);

  /* Vorfaktor */
  factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
  factor *= exp(mp_real(0.25) * pi2 * b * b);
  
  dfactor = 2.0 * (b_in/pow(dpi,1.5)) * exp(b_in*b_in*dpi*dpi*.25);

  cout << "+++ Precision = " << mpipl << " +++" << endl;
  cout << "b  = " << b ;
  cout << "a  = " << a;
  cout << "g  = " << g;
  cout << "h  = " << h;
  cout << "N  = " << N << endl;
  cout << '\v' ;
  cout << "x1 = " << x1;
  cout << "x2 = " << x2;
  cout << "m  = " << m << endl;
  
  cout << "xmax = " << h * N << endl;
  
  cout << "dfaktor = " <<dfactor << endl;
  cout << "Faktor = " <<factor << endl;


  /* Dateien */

  String data = "GAUSSMOD-N" ;
  data = data + N + "-xb" + x2_in + "-beta" + strdup(fix3(g_in)) + ".out"; 
  ofstream dout(data());
  
  String reg = "GAUSSREG-N" ;
  reg = reg + N + "-xb" + x2_in + "-beta" + strdup(fix3(g_in)) + "-b" + strdup(fix3(b_in)) + ".out"; 
  ofstream out(reg());
  
  
  mp_real mp_E2_Gauss, mp_E2_Gausstheo, mp_Reg, mp_Regtheo, mp_ModGauss, step;
  mp_real mp_diffrel, mp_diffrelreg;
  double E2, E2_Gauss, E2_Gausstheo, Reg, Reg_theod, ModGauss, diffrel, x0d;
  double diffrelreg;
  
  step = (x2 -x1)/mp_real(m);


  if(Control == 1)
    {
	cout << "*** Begin: Read Data ***" << endl; 
      for(int j = 0; j <= m; j++)
	{
	  x0 = x1 + j * step; 
	  mp_E2_Gauss = mpe2gauss(x0,g);
	  mp_E2_Gausstheo = mpe2gausstheo(x0,g);
	  mp_real mpdeltarel = abs(mp_E2_Gauss - mp_E2_Gausstheo)/abs(mp_E2_Gausstheo);
	  mpdeltarel *= mp_real(100.0);
	  
	  x0d = dble(x0);
	  E2 = e2(x0d);
	  E2_Gauss        = dble(mp_E2_Gauss);
	  E2_Gausstheo    = dble(mp_E2_Gausstheo);
	  double deltarel = dble(mpdeltarel);
	  cout << x0 << mp_E2_Gausstheo << mp_E2_Gauss << mpdeltarel << endl;;

	  (ostream&) dout << x0d << '\t' << E2 << '\t' << E2_Gausstheo <<'\t' << E2_Gauss << '\t' << deltarel << endl;
      	}
      cout << "*** End: Read Data ***" << endl; 
    }

  cout << "*** Begin: Regularisation ***" << endl;

  for(int jj = 0; jj <= m; jj++)
    {
      x0 = x1 + jj * step;
      s = x0;
      mp_Reg = mptrapez(mpkern,x0,b,g,s,h,N);
      mp_Reg *= factor;
  
      mp_Regtheo  = Regtheo(x0,b,g); //reine Regularisierte
      mp_ModGauss = GaussDichte(x0,g); //Gaussdichte als Modell
      
      
      mp_diffrel = abs(mp_Reg - mp_Regtheo)/abs(mp_Regtheo);
      mp_diffrel *= mp_real(100.0);
      
      mp_diffrelreg  = abs(mp_Reg - mp_ModGauss)/abs(mp_ModGauss);
      mp_diffrelreg *= mp_real(100.0);

      Reg_theod  = dble(mp_Regtheo);
      Reg        = dble(mp_Reg);
      ModGauss   = dble(mp_ModGauss);
      diffrel    = dble(mp_diffrel);
      diffrelreg = dble(mp_diffrelreg);
      x0d        = dble(x0);

      cout << x0 << mp_Reg << mp_Regtheo << mp_ModGauss << mp_diffrel << mp_diffrelreg<< endl;


      (ostream&) out << x0d << '\t' << ModGauss << '\t' << Reg << '\t' << Reg_theod<< '\t' << diffrel <<'\t' << diffrelreg << endl;
      
    }
  cout << "*** End: Regularisation ***" << endl;
  return 0;
}


