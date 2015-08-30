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

/* COS-FFT-CODE */
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
  /*
  double y_,b2_;
  b2_ = b_ * b_;
  y_ = b_/sqrtPi * exp(- b2_ * x_* x_); */
  
  double p,q,y_;
  p = exp(-x_);
  q = 1 + exp(-2.0*x_);
  return y_ = p/q;
}

double F(double x_ , double b_)
{
  //  double y_, b2_;
  //  b2_ = b_ * b_;
  //  y_ = .5 *exp(-.25 * x_ * x_ /b2_);
  double yinv,y_;
  yinv = (2.0/Pi)*cosh(.5*Pi*x_);

  return y_ = 1/(2*yinv);
}

double Kern2(double x)
{
  double y;
  y = (4.0/Pi)*cosh(.5*Pi*x);

  return y;
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
  ofstream out("COS_FFT.out");
  ofstream fout("DATEN.out");
  ofstream finvout("INV_COS_FTT.out");
  out      << setiosflags(ios::uppercase);
  fout     << setiosflags(ios::uppercase);
  finvout  << setiosflags(ios::uppercase);

  cout << setiosflags(ios::uppercase);
  
  
 static unsigned long n;
 double xb, x0, b; // xb = obere Grenze 
 
 cin >> x0 >> b ;
 cin >> xb >> n;
 
 n = pow(2,n);
 //n=2<<18;
 double daten[n+1]; //Datenarray, " C++ build-in type"
 static int nout = n/2;

 // cout << "x0 = " << x0 << endl;
 cout << "Parameter b         = " << b  << endl;
 cout << "Cutoff distanvce xb = " << xb << endl;
 cout << "No. of grides n     = " << n  << endl;
 cout << "max nout  = " << nout  << endl;
 cout << '\v';
 cout << "Grid mesh in the r-space dr = " << xb/n << endl;
 cout << "Grid mesh in the k-space dk = " << Pi/xb << endl;

 cout << '\v';

 /* COS-TRANSFORMATIOMN */

 for(int i = 0; i < n+1;i++) 
   {
   daten[i] = f(xb*i/n,b); 
   fout << i << '\t' << xb*i/n << '\t' << daten[i] << endl;
   }


 cout << "+++ Begin: Cos-Transformation & Write Data +++" << endl;
 dcosft1(daten-1,n);

 for (int i = 0; i < n+1; i++)
   {
     daten[i] *= xb/n; // h = xb/n 
   }
 
 /* Ausgabe der COS-Transformierten und noetige Transformationen */

  /*
 for(int i = 0; i < nout; i++)
   {int j = nstep*(i-nout/2);
   out << j << '\t' << abs(j) << '\t' << j*Pi/xb << '\t' << daten[abs(j)] << endl;
   } */
 
/* nur positive Frequenzen ausgeben */

 double Error_rel;
 double fft_daten;

 for(int i = 0; i < nout; i++)
   { fft_daten = daten[i];
   //     Error_rel = fabs(fft_daten - F(i*Pi/xb,b))/fabs(F(i*Pi/xb,b));
   Error_rel = RelError(F(i*Pi/xb,b),fft_daten);
   out << i << '\t' <<  i*Pi/xb << '\t' << F(i*Pi/xb,b)   << '\t' <<fft_daten << '\t' << Error_rel << endl; 
   }
 

 cout << "*** w_kritisch = " << nout*Pi/xb << '\t' << daten[nout] << '\t' << F(nout*Pi/xb,b)<<" ***" <<  endl;
 
 cout << "+++ End: Cos-Transformation & Write Data +++" << endl;

 cout << '\v';
 cout << "+++ Begin: INV-Cos-Transformation & Write Data +++" << endl;
 
 dcosft1(daten-1,n);

 for (int i = 0; i < n+1; i++)
   {
     daten[i] *= 2.0/xb; 
   }
 
double daten_inv;

 for(int i = 0; i < n; i++)
   { daten_inv = daten[i];
   //     Error_rel = fabs(fft_daten - F(i*Pi/xb,b))/fabs(F(i*Pi/xb,b));
   Error_rel = RelError(f(xb/n*i,b),daten_inv);
   finvout << i << '\t' <<  xb/n*i << '\t' << f(i*xb/n,b)   << '\t' <<daten_inv << '\t' << Error_rel << endl; 
   }

 cout << "+++ End: Cos-Transformation & Write Data +++" << endl;

  return 0;
}


