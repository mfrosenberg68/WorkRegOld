#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "mpmod.H"

mp_real mpf1(mp_real x_)
{ 
  mp_real p;
  mp_real q;
  p =  exp(-x_);
  q = mp_real(1.0) + exp(- 2.0 * x_);
  mp_real y_ = p/q ;
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

mp_real mpIntexakt(mp_real x_, mp_real b_)
{  static mp_real SqrtPi = sqrt(mppic);
   mp_real x2_ = x_ * x_;
   mp_real y = (b_/SqrtPi)*exp(- b_ * b_ * x2_);
return y;
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



int main()
{
  cout << setiosflags(ios::uppercase);
  ofstream out("MP_CONV_SP.out");
  out  << setiosflags(ios::uppercase);


  mp_init();
  mp_real b, x0, s, h, x1, x2;
  double s_in, h_in , bin, x0in, dfactor,dpi, x1_in, x2_in;
  int n,m;
  dpi = 4.0 * atan(1.0);

  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;
  static mp_real factor;
    
  cin >> x0in >> bin;
  cin >> s_in >> h_in >> n;
  cin >> x1_in >> x2_in >> m;

  s = mp_real(s_in);
  h = mp_real(h_in);
  x0 = mp_real(x0in);
  b  = mp_real(bin);
  x1 = mp_real(x1_in);
  x2 = mp_real(x2_in);

  factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
  factor *= exp(mp_real(0.25) * pi2 * b * b);
  dfactor = 2.0 * (bin/pow(dpi,1.5)) * exp(bin*bin*dpi*dpi*.25);

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

  cout << "dfaktor = " <<dfactor << endl;
  cout << "Faktor = " <<factor << endl;

  out << "b = " <<b ;
  out << "x0 = " << x0 << endl;

 
  cout << "Trapez-Integralroutinen" << endl;
  out << "Trapez-Integralroutinen" << endl;

  mp_real mpNIntegral,mpExakt, mpErro_rel;

  mpNIntegral = mptrapez(mpkern,x0,b,s,h,n);
  mpNIntegral *= factor;
  mpExakt = mpIntexakt(x0,b);
  mpErro_rel = abs(mpNIntegral - mpExakt)/abs(mpExakt);
  mpErro_rel *= mp_real(100.0);
 
  cout << "x0 = " << x0;
  out  << "x0 = " << x0;
  cout << mpNIntegral <<mpExakt << mpErro_rel<< endl;
  out << mpNIntegral << mpExakt << mpErro_rel<< endl;

  
  double NIntegral, Exakt, Error_rel;
  NIntegral = dble(mpNIntegral); // kovertiert mp-Zahl nach einer double Zahl
  Exakt = dble(mpExakt);
  Error_rel = dble(mpErro_rel);
  cout <<x0in << '\t' << NIntegral<< '\t' << Exakt <<'\t' << Error_rel << endl ;
  (ostream&)out << x0in << '\t' << NIntegral<< '\t' << Exakt <<'\t' << Error_rel << endl ;

  return 0;
}


