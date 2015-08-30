#include <stdio.h>
#include <iostream.h>  // Input/Output-Operationen
#include <iomanip.h>   // Formtierung Input/Output
#include <fstream.h>   // Einlesen/Schreiben in Files
#include <math.h>

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



double f(double x_, double b_)
{ 
  double y_,b2_;
  b2_ = b_ * b_;
  y_ = b_/sqrtPi * exp(- b2_ * x_ * x_ );
  //  y_ = 1.0/(1.0 + x_*x_);
  return y_; 
  
  /*
  double p,q,y_;
  p = exp(-x_);
  q = 1 + exp(-2.0*x_);
  return y_ = p/q;
  */
}

double FTf(double x, double b)
{

  double y_, b2_;
  b2_ = b * b;
  y_ = exp(-.25 * x * x /b2_);
  //y_ = Pi*exp(-fabs(x));
  return y_;
  /*  
  double yinv,y_;
  yinv = (2.0/Pi)*cosh(.5*Pi*x);

  return y_ = 1/(yinv);
   */
}

int main()
{ 
  ofstream out("FTT_TEST_REAL.out");
  ofstream outfft("FTT_TEST.out");
  ofstream ugout("FTT_TEST_IMAG.out");
  ofstream dout("DATEN-TEST.out");
  ofstream fout("TEST-INV.out");
  ofstream f2out("TEST2-INV.out");


  out  << setiosflags(ios::uppercase);
  outfft  << setiosflags(ios::uppercase);
  ugout << setiosflags(ios::uppercase);
  dout << setiosflags(ios::uppercase);
  fout << setiosflags(ios::uppercase);
  f2out << setiosflags(ios::uppercase);
  cout << setiosflags(ios::uppercase);

  
  
 static unsigned long n;
 double xa, xb, b; // xb = obere Grenze 
 
 
 cin >> b ;
 cin >> xa >> xb >> n;
 
 n = pow(2,n);
 const double dr = (xb-xa)/n;
 const double dk = Pi/xb;
 static unsigned long m = 2*n;
 // const double dk = 1.0/xb;
 // const double a= 1.0/(.5*b);
 
 
 // n=2<<18;

 /* n = 2^m komplexe Datenpunkte im 2*n (!!!) rellem Datenarray  */
 double daten[m+1]; //Datenarray, " C++ build-in type"
 double fftdaten[m+1];
 double x;
 double k;


 cout << "Cutoff distance xa: " << xa << endl;
 cout << "Cutoff distance xb: " << xb << endl;
 // cout << "Parameter a is    : " << a <<endl;
 cout << "Parameter b is    : " << b <<endl;
 cout << "No. of grids      : " << n << endl;
 cout << "Grid mesh in the r-space dr = " << dr << endl;
 cout << "Grid mesh in the k-space dk = " << dk << endl;


 cout << "b  = " << b  << endl;
 cout << "xa = " << xa << endl;
 cout << "xb = " << xb << endl;
 cout << "m  = " << m  << endl;
 cout << "n  = " << n  << endl;
 cout << '\v';

 /* Datenarray auffuellen */

 /* Imaginaerteil der Daten: ungerade Indizes */
 for(unsigned long i = 0; i < (m+1)/2; i++)
   {unsigned long j;
   j = 2*i + 1;
   daten[j] = 0.0;
   }


 /* Realteil der Daten: gerade Indizes */
 
 for(unsigned long i = 0; i < (m+1)/2; i++)
   {unsigned long j;
   j = 2*i;
   x = xa +i*dr;
   daten[j] = f((x+2.0),b); 
   //dout << x << '\t' << daten[j] <<endl;
   }


 /* Daten zwecks Kontrolle ausgeben; Format:*/
 /* x, Realteil, Imaginaerteil */ 

 for(int i = 0; i < (m+1)/2;i++) 
   {
     unsigned long j;
     unsigned long jj;
     j=2*i;
     jj= 2*i+1;
     x = xa +i*dr;
     //x[i] = i*dr;
     dout << x << '\t' << j << '\t' <<  daten[j] << '\t' << jj << '\t' << daten[jj] <<endl;
   }

 /* FTT-Algorhitmus fuer reelle Funktionen aufrufen */

 cout <<"+++ Begin: FFT +++" << endl;
 dfour1(daten-1,n,1); /*"y-1" heißt: Feld beginnt bei 0 und nicht bei 1 num_rec; Pointeruebergabe */ 
 cout <<"+++ Ende: FFT +++" << endl; 

 /* Umskalierung des Ergebnisses */

 for (unsigned long i = 0; i < m+1; i++)
   {
     daten[i] *= dr; // h = xb/n 
     outfft << i << '\t' << daten[i] << endl;
   }

 /* Ausgabe der FFT-Transformierten */

 /* negative Frequnzen */
 double kmax = .5*n*dk;

 cout << "+++ kmax = " << kmax << endl;


 for( unsigned long i = 0; i < (n-1)/2; i++)
   {unsigned long j,jj;
   
   k = (i+1)*dk - kmax;
   j = 2*i + (n+2);
   jj = j+1;
   fftdaten[j] = daten[j]*FTf(k,b);
   fftdaten[jj] = daten[jj]*FTf(k,b);
   // cout << i << '\t' << j << '\t' << jj << endl;
   out << j << '\t' << k << '\t' << daten[j] << '\t' << fftdaten[j] << '\t' << FTf(k,b) << endl ;
   ugout << jj << '\t' << k << '\t' << daten[jj]  << endl;
   }
 


 /* positive Frequenzen */ 

 for(unsigned long i = 0; i <= n/2; i++)
   {unsigned long j,jj;
   k = i*dk;
   j = 2*i;
   jj = 2*i +1;
   fftdaten[j] = daten[j]*FTf(k,b);
   fftdaten[jj] = daten[jj]*FTf(k,b);
   out << j << '\t' << k << '\t' << daten[j] << '\t' << fftdaten[j] << '\t' << FTf(k,b) << endl ;
   ugout << jj << '\t' << k << '\t' << daten[jj]  << endl;
   }


 /* Inversion */
  dfour1(daten-1,n,-1);
  dfour1(fftdaten-1,n,-1);

 for (unsigned long i = 0; i < m+1; i++)  daten[i] *= 1.0/(xb-xa); // h = xb/n 
 for (unsigned long i = 0; i < m+1; i++)  fftdaten[i] *= 1.0/(xb-xa); // h = xb/n 
 
  for(int i = 0; i < (m+1)/2;i++) 
    {
      unsigned long j;
      unsigned long jj;
      j=2*i;
      jj= 2*i+1;
      x = xa +i*dr;
      //x[i] = i*dr;
      fout << x << '\t' <<  daten[j] << '\t' << daten[jj] << '\t' << f(x,b) << '\t' << daten[j] - f(x,b) <<  endl;
      f2out << x << '\t' <<  fftdaten[j] << '\t' << fftdaten[jj] << '\t' << f(x,b) << '\t' << fftdaten[j] - f(x,b) <<  endl;
   }
  
  /*
  for( unsigned long i = 0; i < (n-1)/2; i++)
   {unsigned long j,jj;
   x = xa+ (i+1)*dr ;
   j = 2*i + (n+2);
   jj = j+1;
   f2out << x << '\t' <<  fftdaten[j] << '\t' << fftdaten[jj] << '\t' << f(x,b) << '\t' << fftdaten[j] - f(x,b) <<  endl;
   }
 


  // positive Frequenzen 

 for(unsigned long i = 0; i <= n/2; i++)
   {unsigned long j,jj;
   x = i*dr;
   j = 2*i;
   jj = 2*i +1;
   f2out << x << '\t' <<  fftdaten[j] << '\t' << fftdaten[jj] << '\t' << f(x,b) << '\t' << fftdaten[j] - f(x,b) <<  endl;
   }
   */

  return 0;
}


