#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "nrutil.h"
#include "mpmod.H"     // multiprecision Bibliothek


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
  
  for(int i = 0; i < N; i++)
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
  
  for(int i = 0; i < N; i++)
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


void intlinear(double* x, double* y, int n, double* y2)
{
  double h;
  for(int i = 1; i <= n; i++)
    {
      h = x[i+1] - x[i];
      if (h == 0.0) nrerror("Bad x input to routine intlinear");    
      y2[i] = (y[i+1] - y[i])/h;
    }
}


double linint(double* xa,double* ya, double* y2a, int n, double x)
{
  int klo,khi,k;
  double h,y;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) nrerror("Bad xa input to routine linint");

  y = ya[klo] + y2a[klo] * (x - xa[klo]) ;
  
  return y;
}

double e2data(double x)
{
  double p,q,y,xi;
  //  xi = Xi3double(x);
  p = exp(-x);
  q = 1.0 + exp(- 2.0 * x);
  y = p/q;
  //  y *= (1.0 + 0.1 * xi);
  return y;
}

double e2(double x)
{
  double p,q,y;
  p = exp(-x);
  q = 1.0 + exp(- 2.0 * x);
  y = p/q;
  return y;
}

/* Gaussscher Regularisierungskernel */

mp_real mpregkerne2(mp_real x_, mp_real b_)
{
  static mp_real pi = mppic ;
  mp_real _y = exp(-b_ * b_ * x_ * x_ ) * cos(pi * b_ * b_ * x_ );
  return _y;
}

mp_real mpRegexakt(mp_real x, mp_real b)
{
  static mp_real SqrtPi = sqrt(mppic);
  mp_real y = (b/SqrtPi) * exp(- b * b * x * x);
  return y;
}


int main()
{
  cout  << setiosflags(ios::uppercase);
  
  ofstream out1("Daten_ErrorGauss.out");
  out1 << setiosflags(ios::uppercase);
  
  ofstream out("LIP_Error.out");
  out << setiosflags(ios::uppercase);

  ofstream outs("Y2D_Error.out");
  outs  << setiosflags(ios::uppercase);
  
  ofstream rout("REGLIP_ErrorGauss.out");
  rout  << setiosflags(ios::uppercase);

  ofstream lout("L2-NormGauss.out");
  lout  << setiosflags(ios::uppercase);
  
  cout << "*** Error: Imaginaerteilgleichung ***" << endl;
  cout << "*** Filter: Gauss ***" << endl;

  mp_init();
  cout << "*** Precision = " << mpipl << " ***" << endl;
  

  int M, N, MD;
  double xa, xb, x1, x2, xr1, xr2, h, step, step2, b1, dfactor, dpi;
  mp_real b, mpx1, mpx2, mpstep, mpx;

  dpi = 4.0 * atan(1.0);

  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;
  static mp_real factor;


  cin >> b1;
  cin >> xa >> xb >> MD;
  cin >> x1 >> x2 >> M;  
  cin >> xr1 >> xr2 >> N ;
  
  cout << '\v';

  /* KERNELARRY ANFANG */

  mp_real* const mpregKern = new mp_real [M+1];
  b = mp_real(b1);
  mpx1 = mp_real(x1);
  mpx2 = mp_real(x2);

  mpstep = (mpx2 - mpx1)/mp_real(M);  
  
  /* Vorfaktor */
  factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
  factor *= exp(mp_real(0.25) * pi2 * b * b);
  
  dfactor = 2.0 * (b1/pow(dpi,1.5)) * exp(b1*b1*dpi*dpi*.25);

  cout << "*** Begin: Init. Kernel *** " << endl;
    
  cout << "b1  = " << b1  << endl;
  cout << "b  = " <<  b ;
  cout << "x1  = " << x1 << endl;  /* Anfangs-/Endpunkt der Regularisierung */
  cout << "x2  = " << x2 << endl;
  cout << "M   = " << M  << endl; 
  cout << "h   = " << mpstep;
  cout << "factor = " << factor ;
  cout << "dfactor = " << dfactor << endl;

  for(int i = 0; i <= M; i++) 
    {
      mpx = mpx1 + i * mpstep;
      mpregKern[i] = mpregkerne2(mpx, b);
    }

  cout << "*** End: Init. Kernel *** " << endl;

  /* KERNELARRY ENDE */

  double x0, datlinpol, Error_splint, DataError, e2core;

  /* Index der Vektoren beginnt mit 1 und *nicht* mit 0 */
  
  double* y2d = dvector(1,MD);  /* 2.Ableitung der SPL-Funktion */
  double* data = dvector(1,MD); /* Datenvektor */
  double* x = dvector(1,MD);    /* x-Vektor */
  
  /* Datenarrays fuellen */

  h = (xb - xa)/(MD - 1);
  double rho = MD/(xb -xa);

  cout << '\v';
  cout << "*** Begin: Init. Data-Array ***" << endl;
  cout << "xa  = " << xa << endl;   /* Anfangs-/Endpunkte der Daten */
  cout << "xb  = " << xb << endl;
  cout << "MD  = " << MD  << endl;  /* Zahl der Datenpunkte */
  cout << "h   = " << h  << endl;
  cout << "rho = " << rho  << endl;
  x[1] = xa;
  data[1] = e2data(xa);

  for(int i = 2; i <= MD ; i++)
    {
      x[i] = xa + (i-1) * h;
      data[i] = e2data(x[i]);
    }

  x[MD] = xb;
  data[MD] = e2data(xb);
  
  cout << "*** End: Init. Data-Array ***" << endl;
  cout << '\v' ;

  /* Kontrollausgabe, Datenausgabe */
  for(int i = 1; i <= MD; i++)
    {
      DataError = fabs(data[i] - e2(x[i]))/fabs(e2(x[i]));
      DataError *= 100.0;
      out1 << i << '\t' << x[i] << '\t' << data[i] << '\t' << e2(x[i]) << '\t' << DataError << endl;
    }

  cout << "*** Begin: Interpol. *** " << endl;

  intlinear(x,data,MD,y2d);

  for(int i = 1; i <= MD; i++) outs << i << '\t'  << y2d[i] << endl;

  cout << "*** End:   Interpol. *** " << endl; 
  cout << '\v';

  
  /* L2-Norm des Datenfehlers, incl. Interpolationsfehler */
  double epsx0, eps2data;
  epsx0 = 0.0;

  for(int j = 0; j <= M; j++)
    {
      step = dble(mpstep);
      x0 = x1 + j * step;
      
      datlinpol = linint(x,data,y2d,MD,x0);
      e2core  = e2(x0);
      Error_splint = fabs(datlinpol - e2core)/fabs(e2core);
      Error_splint *= 100.0;

      eps2data = (datlinpol - e2core)*(datlinpol - e2core);      
      epsx0 += eps2data;

      out << x0 << '\t' <<  datlinpol << '\t' << e2core << '\t' << Error_splint << endl;
      
    }
  epsx0 *= step;

  /* REGULARISIERUNG */
 
  double xi;

  mp_real mp_regsum, mpdata, mpError, mpreg;
  double regsum, Reg;
  double Error;

  step2 = (xr2 - xr1)/N;

  cout <<"*** Begin: Error Regularisation ***" << endl;
  cout << "b1  = " << b1  << endl;
  cout << "b  = " <<  b ;
  cout << "x1 = " << xr1 << endl;
  cout << "x2 = " << xr2 << endl;
  cout << "N  = " << N  << endl;
  cout << "h  = " << step2 << endl;
  cout << "*** L2-Datenfehler = " << epsx0 << endl;

  cout << '\v' ;

  mpError = mp_real(0.0); //L2 - Norm des Fehlers

  for(int jj = 0; jj <= N; jj++)
  {
    mp_regsum = mp_real(0.0);
    x0 = xr1 + jj * step2;
    
    for(int j = 0; j <= M; j++)
      {
	xi = x1 + j * step;
	xi = x0 - xi;
	/* Differenz: Interpolierte Daten - exakte Funktion */
	datlinpol = linint(x,data,y2d,MD,xi) - e2(xi);
	mpdata = mp_real(datlinpol);
	mp_regsum += mpregKern[j] * mpdata;
	
      }
    mp_regsum *= mpstep;
    mp_regsum *= factor;

    /* Berechnung des quadrates der L2-Norm des Fehlereinflusses */
    //mpError += step2*exp(x0)*mp_regsum*mp_regsum;

    /* x0 = 0, entspricht nicht! der Norm im urspr. Raum!!! 
     Hier ersteinmal nur zu Testzwecken! */
    mpError += step2*mp_regsum*mp_regsum; 

    mpreg  = mpRegexakt(mp_real(x0),b);
    regsum = dble(mp_regsum);
    Reg    = dble(mpreg);

    cout << x0    << endl;
    cout << mpreg << endl;
    //    cout << mpreg << mp_regsum << mpError << endl;

    Error = dble(mpError);
  
    rout << x0 << '\t' << regsum << '\t' << Reg << '\t' << Reg + regsum <<'\t' << b1 << endl; 
    
  }
  cout << "*** L2-Norm mpError ***" << endl; 
  cout << mpError << endl;
 
  cout << "*** End: Error Regularisation ***" << endl;
  cout << '\v';
  cout << "*** L2-Datenfehler = " << epsx0 << " ***" << endl;
  cout << "*** L2-Norm Error  = " << Error << " ***" << endl; 

  lout << "b " << '\t' << "L2-Datenfehler" << '\t' << " L2-Error" << endl;
  lout << b1 << '\t' << epsx0 << '\t' << Error << endl;
 
  free_dvector(y2d,1,MD);
  free_dvector(data,1,MD);
  free_dvector(x,1,MD);
  delete [] mpregKern;

  return 0;
}
