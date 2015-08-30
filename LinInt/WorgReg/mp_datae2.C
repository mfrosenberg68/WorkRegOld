#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "mpmod.H"
#include "data_mp.H"


int main()
{
  cout  << setiosflags(ios::uppercase);
  ofstream fout("Data.out");
  fout  << setiosflags(ios::uppercase);
  
  mp_init();
  
  double x1_in, x2_in;
  int m;
  
  mp_real b, x0, step, x1, x2;

  cin >> x1_in >> x2_in >> m;

  x1 = mp_real(x1_in);
  x2 = mp_real(x2_in);

  step = (x2 - x1)/m;
 
  cout << "+++ Precision = " << mpipl << " +++" << endl;
  cout << "x1 = " << x1;
  cout << "x2 = " << x2;
  cout << "m  = " << m << endl;
  

  for(int j = 0; j <= m; j++)
    {
      x0 = x1 + j * step;
      mp_real datamp = mpe2(x0);
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

  return 0;
}





