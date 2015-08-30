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
  /*
  ofstream dout("RegE2_MP.out");
  dout  << setiosflags(ios::uppercase);
  */
  mp_init();
  
  double s_in, h_in , b1_in, b2_in, db, dfactor, dpi, x1_in, x2_in;
  int N, Nb, m;

  dpi = 4.0 * atan(1.0);
  
  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;
  static mp_real factor;
  
  mp_real b1, b, b2, x0, s, h, x1, x2, bstep;

  /* Parameter als doubles einlesen */
  cin >> b1_in >> b2_in >> Nb;
  cin >> s_in >> h_in >> N;
  cin >> x1_in >> x2_in >> m;

  /* Konversion: double -> mp_real */
  b1  = mp_real(b1_in);
  b2  = mp_real(b2_in);
  s  = mp_real(s_in);
  h  = mp_real(h_in);
  x1 = mp_real(x1_in);
  x2 = mp_real(x2_in);

  mp_real step = (x2 - x1)/mp_real(m);
  bstep =(b2 -b1)/mp_real(Nb);

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
  cout << "xh  = " << step;
  cout << "xmax = " << h * N << endl;
  
//  cout << "dfaktor = " <<dfactor << endl;
//  cout << "Faktor = " <<factor << endl;
  

  mp_real Regpure_mp, Regtheo_mp, RelNumerror_mp, mpL2_RegError;
//  mp_real RegSol_mp, RelError_mp;
  
  
//  double xod, RegSol, Regtheo, RelError;

  String l2 = "L2-COMLEX-ERROR-MD";
  l2 = l2 + strdup(fix3(N)) + "-xb:" + strdup(fix3(x2_in)) +"-Nb" +strdup(fix3(Nb))+ "-b1:" + strdup(fix3(b1_in)) + "-b2:" + strdup(fix3(b2_in)) + ".out";
  ofstream lout(l2());
  lout << setiosflags(ios::uppercase);

  cout << " *** Beginn PureRegularisation *** " << endl;
  
  for(int jj=0; jj <=Nb; jj++)
  {
      b = b1 + jj * bstep;
      db = dble(b);

      String reg = "GAUSSREG-N" ;
      reg = reg + N + "-xb" + x2_in + "-b" + strdup(fix3(db)) + ".out"; 
      ofstream out(reg());

      /* L2-Norm des Regularsierungsfehlers */
      mpL2_RegError = mp_real(0.0);


        /* Vorfaktor */
      factor = mp_real(2.0)*b/power(pi,mp_real(1.5));
      factor *= exp(mp_real(0.25) * pi2 * b * b);
  
      dfactor = 2.0 * (db/pow(dpi,1.5)) * exp(db*db*dpi*dpi*.25);
      
      cout << " +++ b = " << db << " ;  j = " << jj << " +++" << endl;
      cout << '\v' ;

      for(int j = 0; j <= m; j++)
      {
	 
	 x0 = x1 + j * step;
	 s  = x0; 
	 Regpure_mp  = mptrapez(mpkernpurereg,x0,b,s,h,N);
	 Regpure_mp *= factor;
      
	 Regtheo_mp  = mpRegexakt(x0,b); 
	 
	 RelNumerror_mp = abs(Regpure_mp - Regtheo_mp)/abs(Regtheo_mp);
	 RelNumerror_mp *= 100.0;
	 mpL2_RegError += abs((Regpure_mp - Regtheo_mp)* (Regpure_mp - Regtheo_mp)) ;

	 cout << "x0 = " << x0;
	 cout << Regpure_mp << Regtheo_mp << RelNumerror_mp;
	 cout << '\v' ;
      
	 (ostream&) out << dble(x0) << '\t' << dble(Regtheo_mp) << '\t' << dble(Regpure_mp) << '\t' << dble(RelNumerror_mp) << endl;      
      }
      mpL2_RegError*=step;
      cout << "L2-RegError = " << dble(mpL2_RegError) << endl; 
      cout << '\v';
(ostream&) lout <<  "b = " << db << '\t' <<  " : L2-RegError = " << dble(mpL2_RegError) << endl; 
      
  }
  cout << " *** End PureRegularisation *** " << endl;


  return 0;
}


