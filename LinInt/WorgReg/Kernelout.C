#include <unistd.h>    
#include <iomanip.h>      // Formtierung Input/Output
#include <fstream.h>      // Einlesen/Schreiben in Files
#include <iostream.h>     // Input/Output-Operationen
#include <math.h>         // Default Mathematik-Bibliothek
//#include "nrutil.h"     //NRC Utilities
#include "tnumarray.h"    //DOLIB TemplatArrays :-)
#include "tarray.h"       //DOLIB DoubleArrayx :-)
#include "mpmod.H"        // multiprecision Bibliothek
#include "functionc.H"
#include "stringc.h"
#include "vtype.h"
#include "qc_utilities.H"



int main()
{
  cout  << setiosflags(ios::uppercase);

  int controll_write_kernel = 1;
  /*
  ofstream kout("KERNEL2.out");
  kout << setiosflags(ios::uppercase);
  */

  /* Multiprecisionbibliothek initialisieren */
  mp_init();
  mp_real pi = mppic; // pi aus der Bibliothek !!!
  mp_real pi2 = pi * pi;

  static double  dpi = 4.0 * atan(1.0); // pi als double 
  
  cout << "*** Regularisierende Gauß-Kerne für E1 und E2 ***" << endl;
  cout << "*** Precision = " << mpipl << " ***" << endl;
  
  /* PARAMTER */
  /* Regularisierungsparameter */
  double db1, db2;
  mp_real b1, b2, b;
  int MB;
  static mp_real mpbstep;

  /* Paramter der Regularisierung (Integration)*/
  int M;               // Zahl der Stützstellen

  double x1, x2;       //Anfangs-/Endpunkt
  mp_real mpx1, mpx2;
  static mp_real mpstepreg; //Schrittweite

  /* Parameter einlesen */
  cin >> db1 >> db2 >> MB;
  cin >> x1 >> x2 >> M;  


  /* Konversionen der Parameter: double -> mp_real */
  mpx1 = mp_real(x1);
  mpx2 = mp_real(x2);
  b1    = mp_real(db1);
  b2    = mp_real(db2);
  
  mpbstep =(b2 - b1)/mp_real(MB);

  /* Berechnung weiterer Paramter */
  /* Schrittweite der Integration in der Regularisierung */
  mpstepreg = (mpx2 - mpx1)/mp_real(M);  



  cout << "db1 = " << db1  << endl;
  cout << "b1  = " <<  b1 ;
  cout << "db2 = " << db2  << endl;
  cout << "b2  = " <<  b2 ;
  cout << "x1  = " << x1 << endl;  /* Anfangs-/Endpunkt der Regularisierung */
  cout << "x2  = " << x2 << endl;
  cout << "M   = " << M  << endl; 
  cout << "h   = " << mpstepreg;
  cout << '\v';
  /* ---------------------------------------------------*/

  /* *** Regularisierungskerne und deren Vorfaktoren *** */

  /*  Deklarationen */
  static mp_real factorE1;
  static mp_real factorE2;
  
  NUMARRAY<mp_real> mpregKernE1(M+1);
  NUMARRAY<mp_real> mpregKernE2(M+1);

  //  D_NumArray<mp_real> mpregKernE1(M+1);
  //  D_NumArray<mp_real> mpregKernE2(M+1);

  cout <<"*** Begin: Write Regularisationkernels ***" << endl;    
  cout << '\v';

  for(int j=0; j<=MB; j++)
  {
      b= b1 + j*mpbstep; 

      double db = dble(b);
      String regdatei1 = "KERNE-M";
      regdatei1 = regdatei1 + M + "-b" + strdup(fix3(db)) + ".out";
      ofstream kout(regdatei1());
      kout << setiosflags(ios::uppercase);
  

  /* Definition/Berechnung der Vorfaktoren der Regularisierungskerne */
  factorE1 = mp_real(2.0)*b/power(pi,mp_real(1.5));
  factorE1 *= exp(mp_real(0.25) * pi2 * b * b);
  
  factorE2 = mp_real(2.0)*b/power(pi,mp_real(1.5));
  factorE2 *= exp(mp_real(0.25) * pi2 * b * b);
  

  /* KERNELARRAYs auffuellen */


  cout << " *** b = "  << dble(b) << " ***" << endl;
  cout << "FactorE1  = " << factorE1 ;
  cout << "FactorE2  = " << factorE2 ;
  cout << '\v';
  mp_real mpx;
  
  /* Kontrolle: Regularisierungskern ausgeben */
  if(controll_write_kernel == 1)
    {
	/*
	  ofstream kout("KERNEL2.out");
	  kout << setiosflags(ios::uppercase);

	*/

	//  cout << "*** Kernel Realteil ***" << endl;
	//kout << "*** Kernel Realteil ***" << endl;
      for(int i = 0; i <= M; i++)
	{
	  mpx = mpx1 + i * mpstepreg;
	  //  cout << mpx << mpregkerne1(mpx,b1) << mpregKernE1[i] << endl;
//	  if(controllkernel != mp_real(0.0)) error("Mismatch: Kernel and Kernelarray");  	  
	  kout << dble(mpx) << '\t' <<  dble(factorE1*mpregkerne1(mpx,b)) << '\t' <<  dble(factorE2*mpregkerne2(mpx,b)) << endl;
	}
      
    }
  }
      cout <<"*** End: Write Regularisationkernels ***" << endl;
  
  return 0;
}
