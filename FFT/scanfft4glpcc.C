#include <stdio.h>
#include <iostream.h>  // Input/Output-Operationen
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <math.h>
#include "complex.H"
#include "function.H"
#include "datacc.H"
#include "random.H"
#include "tarray.h"       //DOLIB DoubleArrayx :-)
#include "vtype.h"
#include "qc_utilities.H"


static double Pi = 4.0 * atan(1.0);
static double& pi = Pi;
 
static double sqrtPi = sqrt(4.0 * atan(1.0)) ; // Wurzel von Pi

/* FFT - Algorhitmus */
void swap(double& x,double &y)
{ double z=y; y=x;x=z;}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void dfour1(double *data, unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			swap(data[j],data[i]);
			swap(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP


/* FFT fuer rein reelle Funktionen / Arrays */
void drealft(double *data, unsigned long n, int isign)
{
  //void dfour1();
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		dfour1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		dfour1(data,n>>1,-1);
	}
}

void dcosft1(double *y,int n)
{
	int j,n2;
	double sum,y1,y2;
	double theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;

	theta=pi/n;
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	sum=0.5*(y[1]-y[n+1]);
	y[1]=0.5*(y[1]+y[n+1]);
	n2=n+2;
	for (j=2;j<=(n>>1);j++) {
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=0.5*(y[j]+y[n2-j]);
		y2=(y[j]-y[n2-j]);
		y[j]=y1-wi*y2;
		y[n2-j]=y1+wi*y2;
		sum += wr*y2;
	}
	drealft(y,n,1);
	y[n+1]=y[2];
	y[2]=sum;
	for (j=4;j<=n;j+=2) {
		sum += y[j];
		y[j]=sum;
	}
}


/* FUNKTIONEN */
double F(double x, double b)
{
  double yinv,y_,b2;
  b2 = b*b;
  yinv = (2.0/Pi)*cosh(.5*Pi*x);

  return y_ = exp(-.25 * x * x /b2)/yinv;
}

/*
double TPKern(double x, double b, int p)
{
  double y,kinv,k,b2,M;
  b2 = b*b;
  kinv = (2.0/Pi)*cosh(.5*Pi*x);
  k = 1.0/(kinv);
  M = 1.0+pow(x,2*p);
  return y = k/(k*k + b2*M) ;
}
*/

double Kern2(double x)
{
  double y;
  y = (2.0/Pi)*cosh(.5*Pi*x);

  return y;
}


/*
double ColeDavidson(double x_, double b_, double t0)
{
  double factor = sin(pi*b_)/pi;
  double t = exp(-x_);
  //  double t = x_;
  double y_, kern;
  
  if(t < t0)
    {
      kern = t /(t0 - t);
      //    y_ = factor * power(kern,b_)/t;
      y_ = factor * pow(kern,b_);
    }
  else y_ = double(0.0);

  return y_ ;
}
*/

double StabFkt(double x, double a, double b)
{
  double y;
  double a2 = a*a;
    if (fabs(x) < b)
      {
	return y = exp(-.25 * x * x /a2);
      }
    else
      {
	return y = 0.0;
      }

}

double RelError(double a, double b)
{
  if(a == 0)
    {
      return fabs(b);
    }
  else
    {
      return 100.0*(fabs(a-b)/fabs(a));
    }

}


int main()
{
  cout << setiosflags(ios::uppercase);
  
  static unsigned long n;
  int BN;
  double xa, xb, x0, xif, b, beta, tau0; 
  double b1,b2,a;
  
  cin >> a;
  cin >> b1 >> b2 >> BN;
  cin >> beta >> tau0 >> xif;
  cin >> xa >> xb >> n;
  
  n = pow(2,n);

  const double dx = (xb-xa)/n;
  const double dk = Pi/xb;
  const double bstep = (b2 -b1)/BN;
  double x,k;

  static unsigned long m = 2*n;

  /* n = 2^m komplexe Datenpunkte im 2*n (!!!) rellem Datenarray  */
  double* const daten = new double [m+1]; //Datenarray, " C++ build-in type"
  //  const int nout = n/2;
  //  const int nrout = nout/2;
  
  double kmax = .5*n*dk;
  
  cout << " *** Model: Cole-Cole density *** " << endl;
  cout << "Begin Scan Para. b  = " << b1  << endl;
  cout << "End   Scan Para. b  = " << b2  << endl;
  cout << "No. of scangrides   = " << BN << endl;
  cout << "Scanstep            = " << bstep << endl;
  cout << "Gauss-Parameter  a  = " << a  << endl;
  cout << '\v';

  cout << "Parameter beta      = " << beta  << endl;
  cout << "Parameter tau0      = " << tau0  << endl;
  cout << "Cutoff distanvce xa = " << xa << endl;
  cout << "Cutoff distanvce xb = " << xb << endl;
  cout << "No. of grides n     = " << n  << endl;
  cout << "No. of grides m     = " << m  << endl;
  cout << "kmax                = " << kmax << endl;
  cout << "Noiselevel          = " << xif  << endl;
  cout << '\v';
  cout << "Grid mesh in the r-space dx = " << dx << endl;
  cout << "Grid mesh in the k-space dk = " << dk << endl;
  cout << '\v';

  /* Dateien Anfang */

  ofstream out("CC-DATEN-FFT-REAL.out");
  ofstream ugout("CC-DATEN-FFT-IMA.out");
  ofstream rfout("CC-FFT2REG.out");

  String data = "DATEN-CC-beta";
  data = data + strdup(fix3(beta)) + "-xb" + xb +"-N" + n +".out";
  ofstream fout(data());

  String lscan = "SCAN-BN";
  lscan = lscan + BN + "-a" +strdup(fix3(a)) + "-b2" + strdup(fix3(b2)) + "-b2" + strdup(fix3(b2)) + data; 
  
  ofstream lscanout(lscan());
  

  out      << setiosflags(ios::uppercase);
  fout     << setiosflags(ios::uppercase);
  rfout    << setiosflags(ios::uppercase);
  lscanout    << setiosflags(ios::uppercase);

  /* Dateien Ende */

  /* Datenarray auffuellen und Fehler berechnen*/
  double midrelerror = 0.0;
  double L2_Error = 0.0;
  double s2d = 0.0; 
  double Data;
  int j,jj;


  cout << " +++ Begin: Read Data +++" << endl;

  /* Imaginaerteil der Daten: ungerade Indizes
  for(int i = 0; i < (m+1)/2; i++)
    {
    j = 2*i + 1;
    daten[j] = 0.0;
    }
    */

 /* Realteil der Daten: gerade Indizes */
 
  for(int i = 0; i < (m+1)/2;i++) 
    {
      j = 2*i;
      x = xa + i*dx;
      Data = imag(edata_exakt(x,beta,tau0)); 
      daten[j] = Data * (1.0 + xif*Xi1double(xb*i/n));  
      
      jj = 2*i+1; 
      daten[jj] = 0.0;

      fout << j << '\t' << x << '\t' << Data << '\t' << daten[j] << '\t' << RelError(Data,daten[j]) << '\t' << daten[jj] << endl;  
      midrelerror += RelError(Data,daten[j]);
      L2_Error += (Data - daten[j])*(Data - daten[j]);
    }
  midrelerror /= n;
  L2_Error *= dx;

  for(int i=0; i < (m+1)/2;i++)
    {
      j = 2*i;
      x = xa + i*dx;
      Data = imag(edata_exakt(x,beta,tau0)); 
      s2d += (midrelerror - RelError(Data,daten[j]))*(midrelerror - RelError(Data,daten[j])) ; 
    }

  s2d /= (n-1);

  s2d=sqrt(s2d);
  cout << " +++ End: Read Data +++" << endl;
  cout << '\v';

  double boundmid = UpperBound(midrelerror, 1.0/(2.0*b));
  double boundl2  = UpperBound(L2_Error, 1.0/(2.0*b));
  
  
  /* Datenarray auffuellen und Fehler berechnen Ende*/

  cout << "+++ Begin: FFT-Transformation of Data +++" << endl;
  
  /* FTT-Algorhitmus fuer reelle Funktionen aufrufen:
 "y-1" heißt: Feld beginnt bei 0 und nicht bei 1 num_rec; Pointeruebergabe */ 
  
  dfour1(daten-1,n,1); 
 
 /* Umskalierung des Ergebnisses */

  for (int i = 0; i < m+1; i++) daten[i] *= dx; 
  
  cout << "+++ End: FFT-Transformation of Data +++" << endl;
  cout << '\v';

  cout << "*** Begin: Regularisationscan ***" << endl;
  cout << '\v';

  double Error_rel;
  double fft_daten;
  
  double* const regdaten = new double [m+1]; //Datenarray, " C++ build-in type"
  
  double regsol,daten_inv;
  double RegTh;
  double MidRegErrorRel;
  double L2_RegError;

  for(int l = 0; l <= BN; l++)
    {
      b = b1 +l*bstep;
      cout << " *** b = " << b << "; l = " << l <<" ***" << endl; 

      MidRegErrorRel = 0.0;
      L2_RegError = 0.0;
      
      for (unsigned long i = 0; i < m+1; i++) regdaten[i] = 0.0; 

      String regout = "REGSOLUTION-CC-beta";
      regout = regout + strdup(fix3(beta)) + "-xb" + xb +"-N" + n  +"-a" +strdup(fix3(a)) + "-b" + strdup(fix3(b)) + ".out"; 
      ofstream rout(regout());
      
      String l2 = "L2-ERROR-CC-beta";
      l2 = l2 + strdup(fix3(beta)) + "-xb" + xb +"-N" + n +"-a" +strdup(fix3(a)) + "-b" + strdup(fix3(b)) + ".out";
      ofstream lout(l2());

      rout     << setiosflags(ios::uppercase);
      lout     << setiosflags(ios::uppercase);

  /* negative Frequenzen */
  
  for(int i = 0; i < (n-1)/2; i++)
    {
    k = (i+1)*dk -kmax;
    j = 2*i + (n+2);
    jj = j+1;

    /* Integrand der Regularisierung */
    regdaten[j] = StabFkt(k,a,b)*Kern2(k)*daten[j]; //Realteil der reg. FFT
    regdaten[jj] = StabFkt(k,a,b)*Kern2(k)*daten[jj]; //Imaginaerzeil der reg. FFT

    fft_daten = daten[j];
    Error_rel = RelError(F(k,beta),fabs(daten[j]));
    
    out << j << '\t' << k << '\t' << F(k,beta) << '\t' << fabs(daten[j]) << '\t' << fft_daten << '\t' << Error_rel << endl; 
    
    rfout << j << '\t' << k << '\t' << regdaten[j] << '\t' << StabFkt(k,a,b) << '\t' << daten[j] <<  endl;
    
    ugout << jj << '\t' << k << '\t' << daten[jj] << '\t' << fabs(daten[jj]) << endl;
    }
  
  /* positive Frequenzen */

    for(int i = 0; i < n/2; i++)
    {
      k = i*dk;
      j = 2*i;
      jj =2*i+1;
      
    /* Integrand der Regularisierung */
    regdaten[j] = StabFkt(k,a,b)*Kern2(k)*daten[j]; //Realteil der reg. FFT
    regdaten[jj] = StabFkt(k,a,b)*Kern2(k)*daten[jj]; //Imaginaerzeil der reg. FFT

    fft_daten = daten[j];
    Error_rel = RelError(F(k,beta),fabs(daten[j]));
    
    out << j << '\t' << k << '\t' << F(k,beta) << '\t' << fabs(daten[j]) << '\t' << fft_daten << '\t' << Error_rel << endl; 
    
    rfout << j << '\t' << k << '\t' << regdaten[j] << '\t' << StabFkt(k,a,b) << '\t' << daten[j] <<  endl;
    
    ugout << jj << '\t' << k << '\t' << daten[jj] << '\t' << fabs(daten[jj]) << endl;
    }
  
    dfour1(regdaten-1,n,-1);

  for (int i = 0; i < m+1; i++)
    {
      regdaten[i] *= -1.0/(xb-xa); 
    }
  

  for(int i = 0; i < (m+1)/2; i++)
    { 
      j = 2*i;
      jj = 2*i+1;

      x0 = xa + i *dx;
      regsol = regdaten[j];
      
      //RegTh = e2datagauss_exakt(x0,beta);
      RegTh = ColeCole(-x0,beta,tau0); // theor. Regularisierte
      
      Error_rel = RelError(RegTh,regsol);
      MidRegErrorRel += Error_rel;
      L2_RegError += (RegTh - regsol)*(RegTh - regsol);
      
      rout << -x0 << '\t' << x0 << '\t' << RegTh << '\t' << regsol << '\t' << Error_rel <<'\t' << j <<  endl; 
   }

  MidRegErrorRel = MidRegErrorRel/n;
  L2_RegError *= dx;

  double s2 = 0.0;
  
  for(int i=0; i <= (m+1)/2; i++)
    {
      j = 2*i;
      regsol = regdaten[j];
      x0 = xa + i *dx;
      RegTh = ColeCole(-x0,beta,tau0); // theor. Regularisierte
      
      s2 += (MidRegErrorRel - RelError(RegTh,regsol)) *(MidRegErrorRel - RelError(RegTh,regsol)) ; 
      
    }
  s2 /=(n-1);
  
  s2 = sqrt(s2);


  cout <<" +++ Ave. Data RelError = " << midrelerror << " +/- " << s2d <<" +++" << endl; 
  cout <<" +++ Data L2-Error      = " << L2_Error << " +++" << endl;

  cout <<" +++ Mid.Upper-Bound    = " << boundmid << " +++" << endl;
  cout <<" +++ L2-Upper-Bound     = " << boundl2 << " +++" << endl;
  
  cout <<" +++ Ave. Reg. RelError = " << MidRegErrorRel << " +/- " << s2 << " +++" << endl;
  cout <<" +++ Reg. L2-Error      = " << L2_RegError << " +++" << endl;
  cout << '\v';

  lout <<"Ave. Data RelError = " << midrelerror << " +/- " << s2d << endl;
  lout <<"Data L2-Error      = " << L2_Error << endl;
  lout <<"Mid.Upper-Bound    = " << boundmid << endl;
  lout <<"L2-Upper-Bound     = " << boundl2  << endl;
  lout <<"Ave. Reg. RelError = " << MidRegErrorRel <<" +/- " << s2 << endl;
  lout <<"Reg. L2-Error      = " << L2_RegError << endl;
  lscanout << "b = " << b << ": Reg. L2-Error  = " << L2_RegError << endl;
    }
  
  cout << "*** End: Regularisationscan ***" << endl;
  return 0;
}




