#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "mpmod.H"
#include "data_mp.H"
#include "regkernel_mp.H"
#include "mptrapez.H"
#include "tarray.h"       //DOLIB DoubleArrayx :-)
#include "vtype.h"
#include "qc_utilities.H"

int main()
{
  cout << setiosflags(ios::uppercase);

  //ofstream dout("RegE2_MP.out");
  //dout  << setiosflags(ios::uppercase);

  mp_init();
  
  double s_in, h_in, b1_in, b2_in, db, dfactor, dpi, x1_in, x2_in, dxi;
  int N, m, Nb;

  dpi = 4.0 * atan(1.0);
  
  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;
  static mp_real factor;
  
  mp_real b, b1, b2, bstep, x0, s, h, x1, x2, xi;

  /* Parameter als doubles einlesen */
  cin >> b1_in >> b2_in >> Nb;
  cin >> s_in >> h_in >> N;
  cin >> x1_in >> x2_in >> m >> dxi;

  /* Konversion: double -> mp_real */
  b1  = mp_real(b1_in);
  b2  = mp_real(b2_in);
  s  = mp_real(s_in);
  h  = mp_real(h_in);
  x1 = mp_real(x1_in);
  x2 = mp_real(x2_in);
  xi = mp_real(dxi);

  mp_real step = (x2 - x1)/m;
  bstep = (b2 -b1)/mp_real(Nb);


  cout << "+++ Precision = " << mpipl << " +++" << endl;
  cout << "b1  = " << b1 ;
  cout << "b2  = " << b2 ;
  cout << "hb  = " << bstep;
  cout << "s   = " << s;
  cout << "h   = " << h;
  cout << "N   = " << N << endl;
  cout << '\v' ;
  cout << "x1  = " << x1;
  cout << "x2  = " << x2;
  cout << "m   = " << m << endl;
  cout << "xi  = " << dxi << endl;
  cout << "xh  = " << step;
  cout << "xmax = " << h * N << endl;
  
  //cout << "dfaktor = " <<dfactor << endl;
  //cout << "Faktor = " <<factor << endl;
  
  String data = "Data-N" ;
  data = data + strdup(fix3(N)) + "-xb" +strdup(fix3(h_in*N))  + ".out"; 
  ofstream fout(data());
      

  cout << "*** Begin: Write Data ***" << endl;
  for(int j = 0; j <= m; j++)
    {
       
      x0 = (j * h) - (h*N);
      mp_real datamp = mpe2(x0,xi);
      mp_real fexaktmp = mpe2exakt(x0);
      mp_real Errormp = datamp - fexaktmp;
      mp_real Errorrelmp = 100.0 * (abs(Errormp)/abs(fexaktmp));
      double data = dble(datamp);
      double fexakt = dble(fexaktmp);
      double Error = dble(Errormp);
      double Errorrel = dble(Errorrelmp);
      double xd0 = dble(x0);
      
      (ostream&) fout << xd0 << '\t' << fexakt << '\t' << data << '\t' << Error << '\t' << Errorrel << endl;
    }




  cout << "*** End  : Write Data ***" << endl;

  mp_real Regtheo_mp;
  mp_real RegSol_mp, RelError_mp;

  
  double xod, RegSol, Regtheo, RelError;

  /*
  if(Control == 1) 
    {
      cout << " *** Beginn Testrechnung *** " << endl;
      x0 = x1;
      s  = x0; 
      Regpure_mp  = mptrapez(mpkernpurereg,x0,b,s,h,N);
      Regpure_mp *= factor;
      
      Regtheo_mp  = mpRegexakt(x0,b);

      RelNumerror_mp = abs(Regpure_mp - Regtheo_mp)/abs(Regtheo_mp);
      RelNumerror_mp *= 100.0;
      
      cout << "x0 = " << x0;
      cout << Regpure_mp << Regtheo_mp << RelNumerror_mp;
      
      x0 = mp_real(0.0);
      s  = x0; 
      Regpure_mp  = mptrapez(mpkernpurereg,x0,b,s,h,N);
      Regpure_mp *= factor;
      
      Regtheo_mp  = mpRegexakt(x0,b);

      RelNumerror_mp = abs(Regpure_mp - Regtheo_mp)/abs(Regtheo_mp);
      RelNumerror_mp *= 100.0;
      
      cout << "x0 = " << x0;
      cout << Regpure_mp << Regtheo_mp << RelNumerror_mp;

      x0 = x2;
      s  = x0; 
      Regpure_mp  = mptrapez(mpkernpurereg,x0,b,s,h,N);
      Regpure_mp *= factor;
      
      Regtheo_mp  = mpRegexakt(x0,b);

      RelNumerror_mp = abs(Regpure_mp - Regtheo_mp)/abs(Regtheo_mp);
      RelNumerror_mp *= 100.0;
      
      cout << "x0 = " << x0;
      cout << Regpure_mp << Regtheo_mp << RelNumerror_mp;

      cout << " *** End Testrechnung *** " << endl;

    }
  */
  cout << '\v' ;

  cout << " *** Beginn PureRegularisation *** " << endl;
  
  for(int jj=0; jj <=Nb; jj++)
  {
      b = b1 + jj * bstep;
      db = dble(b);
      
      String reg = "RegGaussE2_MP-N" ;
      reg = reg + strdup(fix3(N)) + "-xb" + strdup(fix3(x2_in)) + "-b" + strdup(fix3(db)) + ".out"; 
      ofstream out(reg());
      
      /* Vorfaktor */
      factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
      factor *= exp(mp_real(0.25) * pi2 * b * b);
      dfactor = 2.0 * (db/pow(dpi,1.5)) * exp(db*db*dpi*dpi*.25);
      
      cout << " +++ b = " << db << " ;  j = " << jj << " +++" << endl;
      cout << '\v' ;
      
  for(int j = 0; j <= m; j++)
    {
      x0 = x1 + j * step;
      s = x0;
      
      RegSol_mp   = mptrapezxi(mpkern,x0,b,xi,s,h,N);
      RegSol_mp  *= factor;

      Regtheo_mp  = mpRegexakt(x0,b);
      
      RelError_mp = abs(RegSol_mp - Regtheo_mp)/abs(Regtheo_mp);
      RelError_mp *= 100.0;

      /* Konversion mp_real -> double */
      xod      = dble(x0);
      RegSol   = dble(RegSol_mp);
      Regtheo  = dble(Regtheo_mp);
      RelError = dble(RelError_mp);
      
      cout <<  x0 << Regtheo_mp << RegSol_mp << RelError_mp << endl;
(ostream&) out << xod << '\t' << Regtheo << '\t' << RegSol << '\t' << RelError << endl;
    }

  }

  return 0;
}


