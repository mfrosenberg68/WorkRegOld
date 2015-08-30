#include <stdio.h>
#include <iostream.h>  // Input/Output-Operationen
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <math.h>
#include "function.H"
#include "datac.H"
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

double TPKern(double x, double b, int p)
{
  double y,kinv,k,b2,M;
  b2 = b*b;
  kinv = (2.0/Pi)*cosh(.5*Pi*x);
  k = 1.0/(kinv);
  M = 1.0+pow(x,2*p);
  return y = k/(k*k + b2*M) ;
}


double StabFkt(double x, double b)
{
  double y, b2;
  b2 = b * b;
  return y = exp(-.25 * x * x /b2);
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
  int p;
  double xa, xb, x0, xif, b;
  double beta1, beta2, t1, t2, a; 

  cin >> b >> p;
  cin >> beta1 >> beta2;
  cin >> t1 >> t2;
  cin >> a >> xif;
  cin >> xa >> xb >> n;
  
  n = pow(2,n);
  
  double daten[n+1]; //Datenarray, " C++ build-in type"
  //  const int nout = n/2;
  //  const int nrout = nout/2;
  
  const double dx = (xb-xa)/n;
  const double dk = Pi/xb;
  double x,k;


  cout << " *** Model: Gauss-density *** " << endl;
  cout << "Parameter b         = " << b  << endl;
  cout << "Parameter p         = " << p  << endl;
  cout << "Parameter beta1     = " << beta1  << endl;
  cout << "Parameter beta2     = " << beta2  << endl;
  cout << "Parameter t1        = " << t1  << endl;
  cout << "Parameter t2        = " << t2  << endl;
  cout << "Parameter a1        = " << a  << endl;
  cout << "Parameter a2        = " << 1.0 - a  << endl;
  
  cout << "Cutoff distanvce xa = " << xa << endl;
  cout << "Cutoff distanvce xb = " << xb << endl;
  cout << "No. of grides n     = " << n  << endl;
  //  cout << "max nout            = " << nout  << endl;
  //  cout << "max #reg-out        = " << nrout  << endl;
  cout << "Noiselevel          = " << xif  << endl;
  cout << '\v';
  cout << "Grid mesh in the r-space dx = " << dx << endl;
  cout << "Grid mesh in the k-space dk = " << dk << endl;
  //  cout << "Grid mesh in the rout-space = " << xb/nout << endl;
  cout << '\v';

  /* Dateien Anfang */

  ofstream out("TP-GAUSS-DATEN-FFT-REAL.out");
  ofstream ugout("TP-GAUSS-DATEN-FFT-IMA.out");
  ofstream rfout("GAUSS-FFT2REG.out");

  String data = "TP-DATEN-GAUSS-beta1";
  data = data + strdup(fix3(beta1))+"-"+strdup(fix3(t1)) + "-beta2" + strdup(fix3(beta2))+"-"+strdup(fix3(t2)) + "-xb" + xb +"-N" + n +".out";
  ofstream fout(data());

  String regout = "TP-REGSOLUTION-GAUSS-beta1";
  regout = regout +strdup(fix3(beta1))+"-"+strdup(fix3(t1)) + "-beta2" + strdup(fix3(beta2))+"-"+strdup(fix3(t2)) + "-xb" + xb +"-N" + n  + "-p" + p + "-b" +  strdup(sci3(b)) + ".out"; 
  ofstream rout(regout());
  
  String l2 = "TP-L2-ERROR-GAUSS-beta1";
  l2 = l2 + strdup(fix3(beta1))+"-"+strdup(fix3(t1)) + "-beta2" + strdup(fix3(beta2)) +"-"+strdup(fix3(t2)) + "-xb" + xb +"-N" + n + "-p" + p + "-b" +  strdup(sci3(b)) + ".out";
  ofstream lout(l2());


  out      << setiosflags(ios::uppercase);
  fout     << setiosflags(ios::uppercase);
  rout     << setiosflags(ios::uppercase);
  rfout    << setiosflags(ios::uppercase);

  /* Dateien Ende */

  /* Datenarray auffuellen und Fehler berechnen*/
  double midrelerror = 0.0;
  double L2_Error = 0.0;
  double s2d = 0.0; 
  double Data;

  cout << " +++ Begin: Read Data +++" << endl;
  for(int i = 0; i < n+1;i++) 
    {
      x = xa + i*dx;
      Data = e2data_multi_gauss(x,beta1,t1,beta2,t2,a);
      daten[i] = Data * (1.0 + xif*Xi1double(xb*i/n));  
      
      fout << i << '\t' << x << '\t' << Data << '\t' << daten[i] << '\t' << RelError(Data,daten[i]) << endl;  
      midrelerror += RelError(Data,daten[i]);
      L2_Error += (Data - daten[i])*(Data - daten[i]);
    }
  midrelerror /= n;
  L2_Error *= dx;

  for(int i=0; i < n+1;i++)
    {
      x = xa + i*dx;
      //      Data = e2datagauss_exakt(x,beta);
      Data = e2data_multi_gauss(x,beta1,t1,beta2,t2,a);
      s2d += (midrelerror - RelError(Data,daten[i]))*(midrelerror - RelError(Data,daten[i])) ; 
    }

  s2d /= (n-1);

  s2d=sqrt(s2d);
  cout << " +++ End: Read Data +++" << endl;
  cout << '\v';
  
  /* Datenarray auffuellen und Fehler berechnen Ende*/

  cout << "+++ Begin: FFT-Transformation of Data +++" << endl;
  
  /* FTT-Algorhitmus fuer reelle Funktionen aufrufen:
 "y-1" heißt: Feld beginnt bei 0 und nicht bei 1 num_rec; Pointeruebergabe */ 
  
  drealft(daten-1,n,1); 
 
 /* Umskalierung des Ergebnisses */

  for (unsigned long i = 0; i < n+1; i++) daten[i] *= dx; 

  cout << "*** kmax = " <<.5*(n+1)*dk <<'\t'<< F(.5*(n+1)*dk,beta1)  << '\t' << daten[n] << '\t' << RelError(F(.5*(n+1)*dk,beta1),daten[n])<<" ***" << endl; 

  cout << "+++ End: FFT-Transformation of Data +++" << endl;
  cout << '\v';

  cout << "+++ Begin: Transformation & Write Data +++" << endl;

/* nur positive Frequenzen ausgeben */

  double Error_rel;
  double fft_daten;
  
  double regdaten[n+1]; //Datenarray, " C++ build-in type"
 
  /* Imaginaeteil der FFT = ungerader Index */
  
  for(unsigned long i = 0; i < (n+1)/2; i++)
    {unsigned long j;
    k = i*dk;
    j = 2*i + 1;
    //    regdaten[j] = daten[j];
    regdaten[j] = 0.0;
    
    ugout << j << '\t' << k << '\t' << daten[j] << '\t' << fabs(daten[j]) << endl;
    }
  
  
  /* Realteil der FFT = gerader Index incl. i=0 */
  
  out << 0 << '\t' << 0 << '\t' << F(0,beta1) << '\t' << daten[0] << '\t' << fabs(daten[0]) << '\t' << RelError(F(0,beta1),fabs(daten[0])) << endl; 

  out << 1 << '\t' << .5*(n+1)*dk << '\t' << F(.5*(n+1)*dk,beta1) << '\t' << daten[1] << '\t' << fabs(daten[1]) << '\t' << RelError(F(.5*(n+1)*dk,beta1),fabs(daten[1])) << endl; 

  
  regdaten[0] = TPKern(0.0,b,p)*daten[0] ;
  regdaten[1] = TPKern(.5*(n+1)*dk,b,p) * daten[1];
    
  for(unsigned long i = 1; i <=(n+1)/2; i++)
    { 
      unsigned long j;
      j = 2*i;
      k = i*dk;
      fft_daten = daten[j];
      
      //Integrand der Regularisierung 
      regdaten[j] = TPKern(k,b,p)*daten[j];
      
      Error_rel = RelError(F(k,beta1),fabs(daten[j]));
      out << j << '\t' << k << '\t' << F(k,beta1) << '\t' << fabs(daten[j]) << '\t' << fft_daten << '\t' << Error_rel << endl; 
      
      rfout << j << '\t' << k << '\t' << regdaten[j] << '\t' << TPKern(k,b,p) << '\t' << daten[j] <<  endl;
      // rfout << i << '\t' << k << '\t' << StabFkt(k,b) <<  endl;
    }
  
  cout << "+++ End: Transformation & Write Data +++" << endl;
  cout << '\v';

  cout << "+++ Begin: INV-FFT & Regularisation +++" << endl;
 
  drealft(regdaten-1,n,-1);
  drealft(daten-1,n,-1);

  
  for (int i = 0; i < n+1; i++)
    {
      daten[i] *= 2.0/(xb-xa); 
    }
  

  for (int i = 0; i < n+1; i++)
    {
      regdaten[i] *= 2.0/(xb-xa); 
    }
  
  double regsol,daten_inv;
  double RegTh;
  double MidRegErrorRel = 0.0;
  double L2_RegError    = 0.0;


  //  double dxout = (xb-xa)/nout;

  for(int i = 0; i <= n; i++)
    { 
      regsol = regdaten[i];
      x0 = xa + i *dx;
     
      //RegTh = e2datagauss_exakt(x0,beta);
      //      RegTh = GAUSS(x0,beta); // theor. Regularisierte
      RegTh = a*GAUSS((t1 - x0),beta1) + (1.0-a)* GAUSS((t2 - x0),beta2);
      Error_rel = RelError(RegTh,regsol);
      MidRegErrorRel += Error_rel;
      L2_RegError += (RegTh - regsol)*(RegTh - regsol);
      
      rout << i << '\t' << x0 << '\t' << RegTh << '\t' << regsol << '\t' << Error_rel << endl; 
   }

  MidRegErrorRel = MidRegErrorRel/n;
  L2_RegError *= dx;

  double s2 = 0.0;
  
  for(int i=0; i <= n; i++)
    {
      regsol = regdaten[i];
      x0 = xa + i *dx;
      //      RegTh = GAUSS(x0,beta); // theor. Regularisierte
      RegTh = a*GAUSS((t1 - x0),beta1) + (1.0-a)* GAUSS((t2 - x0),beta2);
      s2 += (MidRegErrorRel - RelError(RegTh,regsol)) *(MidRegErrorRel - RelError(RegTh,regsol)) ; 
      
    }
  s2 /=(n-1);
  
  s2 = sqrt(s2);

  cout << "+++ End: INF-FFT & Regularisation +++" << endl;
  cout << '\v';

  cout <<" +++ Ave. Data RelError = " << midrelerror << " +/- " << s2d <<" +++" << endl; 
  cout <<" +++ Data L2-Error      = " << L2_Error << " +++" << endl;
  
  cout <<" +++ Ave. Reg. RelError = " << MidRegErrorRel << " +/- " << s2 << " +++" << endl;
  cout <<" +++ Reg. L2-Error      = " << L2_RegError << " +++" << endl;


  lout <<"Ave. Data RelError = " << midrelerror << " +/- " << s2d << endl;
  lout <<"Data L2-Error      = " << L2_Error << endl;

  lout <<"Ave. Reg. RelError = " << MidRegErrorRel <<" +/- " << s2 << endl;
  lout <<"Reg. L2-Error      = " << L2_RegError << endl;
  
  
  return 0;
}




