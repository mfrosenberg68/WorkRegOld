#include <unistd.h>    
#include <iomanip.h>      // Formtierung Input/Output
#include <fstream.h>      // Einlesen/Schreiben in Files
#include <iostream.h>     // Input/Output-Operationen
#include <math.h>         // Default Mathematik-Bibliothek
#include "mpmod.H"        // multiprecision Bibliothek

//const DComplex I     = DComplex(0.0,1.0);

mp_complex mptrapez(mp_complex function(mp_real _x, mp_real _x0, mp_real _a), mp_real _x0, mp_real _a, mp_real s_, mp_real h_, int n_) 
{

  mp_complex sum;
  mp_real zu, zd;
  //  cout <<"from trapez: s = " << s_ ;
  sum = function(s_,_x0,_a);
  for (int i = 1; i <= n_ ; i++) {
    zu = s_ + i * h_ ;
    zd = s_ - i * h_ ;
    //    cout <<"from trapez: zu = " << zu;

    sum += function(zu,_x0,_a) + function(zd,_x0,_a);
  }
  
  return sum *= h_ ;
}

mp_complex mpregkernc(mp_real x, mp_real x0, mp_real b)
{
  static mp_real pi   = mppic ;
  static mp_real pi2  = mppic*mppic;
  static mp_complex I = mp_complex(0.0,1.0);
  mp_real b2 = b*b;
  
  mp_complex y = exp(-b2*(x0-x)*(x0-x))*(cos(pi*b2*(x0-x))+I*sin(pi*b2*(x0-x)))*(I*cos(2.0*pi*b2*(x0-x))*sinh(pi2*b2)-sin(2.0*pi*b2*(x0-x))*cosh(pi2*b2));

  return y;

}

mp_real Regtheo(mp_real x_, mp_real b_)
{
  static mp_real Sqrtpi = sqrt(mppic);
  mp_real y_ = b_/Sqrtpi * exp(- b_ * b_ * x_ * x_);
  return y_ ;
}

mp_complex E_data(mp_real x)
{
  const mp_complex I = mp_complex(0.0,1.0);

  mp_complex  g = 1.0/(1.0 - I * exp(-x));
  return g;
}

mp_complex RegIntegrand(mp_real x, mp_real x0, mp_real b)
{
  mp_complex Integrand = mpregkernc(x,x0,b)*E_data(x);
  return Integrand;
}

int main()
{
  mp_init();

  //  const mp_complex I = mp_complex(0.0,1.0);

  static mp_real pi  = mppic;
  static mp_real pi2 = mppic * mppic;

  static mp_real Sqrtpi = sqrt(mppic);
  
  cout  << setiosflags(ios::uppercase);

  ofstream out("REGSOLC.out");
  out << setiosflags(ios::uppercase);

  ofstream kout("KERNELCMP.out");
  kout << setiosflags(ios::uppercase);
  
  ofstream dout("DATACMP.out");
  dout << setiosflags(ios::uppercase);
  
  int M, N;
  double bd;
  double x1d, x2d, xad, xbd;
  
  cin >> bd;
  cin >> xad >> xbd >> N;
  cin >> x1d >> x2d >> M;

  cout << "*** Precision = " << mpipl << " ***" << endl;
  cout << "b  = " << bd << endl;
  cout << "xa = " << xad << endl;
  cout << "xb = " << xbd << endl;
  cout << "x1 = " << x1d << endl;
  cout << "x2 = " << x2d << endl;
  cout << "N  = " << N << endl;
  cout << "M  = " << M << endl;

  mp_real xa   = mp_real(xad);
  mp_real xb   = mp_real(xbd);
  mp_real x1   = mp_real(x1d);
  mp_real x2   = mp_real(x2d);
  mp_real  b   = mp_real(bd);

  mp_real factor;
  factor = b/power(pi,mp_real(1.5))*exp(mp_real(1.25)*pi2*b*b);


  mp_real h    = (xb - xa)/mp_real(N);
  cout << "h  = " << h ;

  mp_real step = (x2 - x1)/mp_real(M);
  cout << "h2 = " << step;

  cout << "Faktor = " << factor << endl;


  mp_real mpx0; 
  //  mp_real mpx, mpx0; 
  //  mp_complex regkern, data;
  //  double Realregkern, Imaregkern;
  //  double RealData, ImaData;
  /*
  cout << " *** write Kernel and Data *** " << endl;
  
  for(int i = 0; i <= N; i++)
    {
      mpx0 = mp_real(0.0);
      mpx = xa + i * h;
      regkern = mpregkernc(mpx,mpx0,b);
      //      cout << mpx << regkern << endl;

      Realregkern = dble(MP_REAL(regkern));
      Imaregkern  = dble(aimag(regkern));

      //      data = E_data(mpx);
      //      cout << mpx << data << endl;

      //      RealData = dble(MP_REAL(data));
      //      ImaData  = dble(aimag(data));

      kout << dble(mpx) << '\t' << Realregkern << '\t' << Imaregkern << '\t' << bd<< endl;

      //      dout << dble(mpx) << '\t' << RealData << '\t' << ImaData << endl;
    }
    */
  mp_complex mpregsum;
  mp_real mpReg, RelError;
  double RealReg, ImaReg, Error;

  cout << " *** Begin: Regularisation *** " << endl;

  for(int j = 0; j <= M; j++)
    {
      mpx0 = x1 + j * step;

      mpregsum = mptrapez(RegIntegrand,mpx0,b,mpx0,.5*h,N);
      mpregsum *=factor;

      mpReg = Regtheo(mpx0,b);
      RelError = abs(mpregsum - mpReg)/abs(mpReg);
      RelError *= mp_real(100.0);

      cout << mpx0 << mpregsum  << endl;
      cout << mpReg << RelError << endl;
      
      RealReg = dble(MP_REAL(mpregsum));
      ImaReg  = dble(aimag(mpregsum));
      Error   = dble(RelError);

      out << dble(mpx0) << '\t' << dble(mpReg) << '\t' << RealReg << '\t' << ImaReg << '\t' << Error << '\t' << bd << endl;
    }
  cout << " *** End: Regularisation *** " << endl;
  
  return 0;
}





