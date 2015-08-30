#include <unistd.h>    
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <iostream.h>  // Input/Output-Operationen
#include <math.h>      // Default Mathematik-Bibliothek
#include "nrutil.h"
#include "mpmod.H"     // multiprecision Bibliothek


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

  cout << "klo = " << klo << endl;
  cout << "khi = " << khi << endl;
  cout << "x["<<klo<<"] = " << xa[klo] << endl;
  cout << "x["<<khi<<"] = " << xa[khi] << endl;
  cout << "x    = " << x << endl; 
  cout << '\v';

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

int main()
{
  cout  << setiosflags(ios::uppercase);
  
  ofstream out1("Daten.out");
  out1 << setiosflags(ios::uppercase);
  
  ofstream rout("Int.out");
  rout << setiosflags(ios::uppercase);


  ofstream outs("Y2D.out");
  outs  << setiosflags(ios::uppercase);
  
  double xa, xb, x1, x2, h, step;
  int MD, RP;

  cin >> xa >> xb >> MD;
  cin >> x1 >> x2 >> RP;
  
   /* Index der Vektoren beginnt mit 1 und *nicht* mit 0 */
  
  double* y2d = dvector(1,MD);  /* 2.Ableitung der SPL-Funktion */
  double* data = dvector(1,MD); /* Datenvektor */
  double* x = dvector(1,MD);    /* x-Vektor */
  
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
//cout << i << '\t' << x[i] << '\t' << data[i] << '\t' << e2data(x[i]) << endl; 
    }

  x[MD] = xb;
  data[MD] = e2data(xb);
  
  cout << "*** End: Init. Data-Array ***" << endl;
  cout << '\v' ;

  /* Kontrollausgabe, Datenausgabe */
  for(int i = 1; i <= MD; i++)
    {
      // DataError = fabs(data[i] - e2(x[i]))/fabs(e2(x[i]));
      //      DataError *= 100.0;
      //out1 << i << '\t' << x[i] << '\t' << data[i] << '\t' << e2(x[i]) << '\t' << DataError << endl;
      out1 << i << '\t' << x[i] << '\t' << data[i] << endl;
    }

  cout << "*** Begin: Interpol. *** " << endl;

  intlinear(x,data,MD,y2d);

  for(int i = 1; i <= MD; i++) outs << i << '\t'  << y2d[i] << endl;

  cout << "*** End:   Interpol. *** " << endl; 
  cout << '\v';

  double xi, datlinpol;

  step = (x2 - x1)/RP;

  for(int i = 0; i <= RP; i++)
    {
      xi = x1 + i * step;
      datlinpol = linint(x,data,y2d,MD,xi);
      rout << i << '\t' << xi << '\t' << datlinpol << endl;
    }

  free_dvector(y2d,1,MD);
  free_dvector(data,1,MD);
  free_dvector(x,1,MD);


  return 0;
}
