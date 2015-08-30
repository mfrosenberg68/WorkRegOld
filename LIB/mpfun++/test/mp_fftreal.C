#include <iostream.h>
#include "mpmod.H"

extern "C" {
    #include <time.h>
    #include <stdlib.h>
    #include <math.h>
}

extern void rcfftz(const int, const int, mp_real*, mp_real*, mp_real*);
extern void crfftz(const int, const int, mp_real*, mp_real*, mp_real*);
extern void cfftz(const int, const int, mp_real*, mp_real*, mp_real*);
extern void fftz1(const int, const int, const int, const mp_real*,
		  mp_real*, mp_real*);
extern void fftz2(const int, const int, const int, const mp_real*,
		  mp_real*, mp_real*);
extern void trans(const int, const int, const mp_real*, mp_real*);

int main()

/*
!   Performs real-to-complex and complex-to-real FFT routines using routines
!   RCFFTZ, CRFFTZ and CFFTZ.  These routines use double precision data, not
!   complex.  Converted to Fortran-90 based multiprecision.

!   David H. Bailey     May 4, 1994
!
!   Translated to C++.
!
!   Siddhartha Chatterjee (RIACS).
!   25 August 1994.
*/

{
    const int m = 10;
    //const int n = (int) pow(2, m);
    const int n = 1024;
    const int it = 10;
    mp_real u[2*n], x[n+2], y[n+2], se;
    double s[n];
    time_t t0, t1;

    //!   Initialize.  The SECOND function is assumed to be the timing function
    //!   on the given computer.  Change here and below if another name is used.

    mp_init();
    double rn = 1.0/n;
    cfftz(0, m, u, x, y);
    time(&t0);

    //!  Perform IT iterations.
    for (int ii = 1; ii <= it; ii++) {
	for (int i = 0; i < n; i++) {
	    s[i] = dble(mp_rand());
	    x[i] = s[i];
	}

	rcfftz(-1, m, u, x, y);
	crfftz(1, m, u, x, y);
	se = 0.0;

	for (i = 0; i < n; i++)
	  se += power(rn*x[i]-s[i], 2);

	se = sqrt(se/n);
	cout << "Iteration " << ii << "   Error = \n";
	cout << se;
	if (se > mpeps) {
	    cout << "*** Error detected\n";
	    exit(1);
	}
    }

    time(&t1);
    double tm = difftime(t1, t0);
    double rt = it*n*5.0*(m+3.0)/tm*1e-6;
    cout << "Test passed.\n";
//    cout << "CPU time = " << tm << '\n';
//    cout << "MFLOPS rate = " << rt << '\n';
    return 0;
}

void rcfftz(const int is, const int m, mp_real* u, mp_real* x,
	    mp_real* y)

/*
!   Performs an N-point real-to-complex FFT, where N = 2^M.  X is both the
!   input and the output data array, and Y is a scratch array.  N real values
!   are input in X, and N/2 + 1 complex values are output in X, with real and
!   imaginary parts separated by N/2 + 1 locations.  Before calling RCFFTZ,
!   the U array must be initialized by calling CFFTZ with IS = 0.  A call to
!   RCFFTZ with IS = 1 (or -1) indicates a call to perform a real-to-complex
!   FFT with positive (or negative) exponentials.  The arrays X and Y must be
!   dimensioned with N/2 + 1 complex or N + 2 real cells.  U must be
!   dimensioned the same as in CFFTZ.  The output values from RCFFTZ are
!   twice as large as the results of a complex-to-complex transform on real
!   data.  M must be at least three.

!   This subroutine employs a technique that converts a real-to-complex FFT
!   to a complex-to-complex FFT.  See Brigham, "The Fast Fourier Transform",
!   p. 169, and Hockney & Jesshope, "Parallel Computers", p. 303, although
!   neither reference is complete and error-free.

!   David H. Bailey     May 4, 1988
*/

{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    //!   Set initial parameters.
    int n = (int) pow(2, m);
    int mx = INT(u[1-1]);
    int nu = (int) pow(2, mx);
    int n2 = n/2;
    int n21 = n2+1;
    int n4 = n/4;
    int ku = n/2;
    int kn = ku+nu;

    //!   Check if input parameters are invalid.
    if ((is != 1 && is != -1) || m < 3 || m > mx) {
	cout << "RCFFTZ: Either U has not been initialized, or else\n";
	cout << "one of the input parameters is invalid ";
	cout << is << ' ' << m << ' ' << mx << '\n';
	exit(1);
    }

    //!   Copy X to Y such that Y(k) = X(2k-1) + i X(2k).
    //!   This loop is vectorizable.
    for (int k = 1; k <= n2; k++) {
	y[k-1] = x[2*k-1-1];
	y[k+n2-1] = x[2*k-1];
    }

    //!   Perform a normal N/2-point FFT on Y.
    cfftz(is, m-1, u, y, x);

    //!   Reconstruct the FFT of X.
    x[1-1] = 2.0*(y[1-1]+y[n21-1]);
    x[n21+1-1] = 0.0;
    x[n4+1-1] = 2.0*y[n4+1-1];
    x[n4+1+n21-1] = 2.0*is*y[n4+n2+1-1];
    x[n21-1] = 2.0*(y[1-1]-y[n21-1]);
    x[n+2-1] = 0.0;

    //!   This loop is vectorizable.
    for (k = 2; k <= n4; k++) {
	mp_real y11 = y[k-1];
	mp_real y12 = y[k+n2-1];
	mp_real y21 = y[n2+2-k-1];
	mp_real y22 = y[n+2-k-1];
	mp_real a1 = y11+y21;
	mp_real a2 = y11-y21;
	mp_real b1 = y12+y22;
	mp_real b2 = y12-y22;
	mp_real u1 = u[k+ku-1];
	mp_real u2 = is*u[k+kn-1];
	mp_real t1 = u1*b1+u2*a2;
	mp_real t2 = -u1*a2+u2*b1;
	x[k-1] = a1+t1;
	x[k+n21-1] = b2+t2;
	x[n2+2-k-1] = a1-t1;
	x[n+3-k-1] = -b2+t2;
    }
}

void crfftz(const int is, const int m, mp_real* u, mp_real* x,
	    mp_real* y)

/*
!   Performs an N-point complex-to-real FFT, where N = 2^M.  X is both the
!   input and the output data array, and Y is a scratch array.  N/2 + 1
!   complex values are input in X, with real and imaginary parts separated by
!   N/2 + 1 locations, and N real values are output.  Before calling CRFFTZ,
!   the U array must be initialized by calling CFFTZ with IS = 0.  A call to
!   CRFFTZ with IS = 1 (or -1) indicates a call to perform a complex-to-real
!   FFT with positive (or negative) exponentials.  The arrays X and Y must be
!   dimensioned with N/2 + 1 complex or N + 2 real cells.  U must be
!   dimensioned the same as in CFFTZ.  M must be at least three.

!   This subroutine employs a technique that converts a complex-to-real FFT
!   to a complex-to-complex FFT.  See Brigham, "The Fast Fourier Transform",
!   p. 169, and Hockney & Jesshope, "Parallel Computers", p. 303, although
!   neither reference is complete and error-free.

!   David H. Bailey     May 4, 1988
*/

{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    //!   Set initial parameters.
    int n = (int) pow(2, m);
    int mx = INT(u[1-1]);
    int nu = (int) pow(2, mx);
    int n2 = n/2;
    int n21 = n2+1;
    int n4 = n/4;
    int ku = n/2;
    int kn = ku+nu;

    //!   Check if input parameters are invalid.
    if ((is != 1 && is != -1) || m < 3 || m > mx) {
	cout << "CRFFTZ: Either U has not been initialized, or else\n";
	cout << "one of the input parameters is invalid ";
	cout << is << ' ' << m << ' ' << mx << '\n';
	exit(1);
    }

    //!   Construct the input to CFFTZ.
    y[1-1] = 0.5*(x[1-1]+x[n21-1]);
    y[n2+1-1] = 0.5*(x[1-1]-x[n21-1]);
    y[n4+1-1] = x[n4+1-1];
    y[n4+n2+1-1] = -is*x[n4+n2+2-1];

    //!      
    //!   This loop is vectorizable.
    for (int k = 2; k <= n4; k++) {
	mp_real x11 = x[k-1];
	mp_real x12 = x[k+n21-1];
	mp_real x21 = x[n2+2-k-1];
	mp_real x22 = x[n+3-k-1];
	mp_real a1 = x11+x21;
	mp_real a2 = x11-x21;
	mp_real b1 = x12+x22;
	mp_real b2 = x12-x22;
	mp_real u1 = u[k+ku-1];
	mp_real u2 = is*u[k+kn-1];
	mp_real t1 = u1*b1+u2*a2;
	mp_real t2 = u1*a2-u2*b1;
	y[k-1] = 0.5*(a1-t1);
	y[k+n2-1] = 0.5*(b2+t2);
	y[n2+2-k-1] = 0.5*(a1+t1);
	y[n+2-k-1] = 0.5*(-b2+t2);
    }

    //!   Perform a normal N/2-point FFT on Y.
    cfftz(is, m-1, u, y, x);

    //!   Copy Y to X such that Y(k) = X(2k-1) + i X(2k).
    //!   This loop is vectorizable.

    for (k = 1; k <= n2; k++) {
	x[2*k-1-1] = y[k-1];
	x[2*k-1] = y[k+n2-1];
    }
}

int min(const int m, const int n) {return m < n ? m : n;}

void cfftz(const int is, const int m, mp_real* u, mp_real* x, mp_real* y)

/*
!   Computes the 2^M-point complex-to-complex FFT of X using an algorithm due
!   to Swarztrauber, coupled with some fast methods for performing power-of-
!   two matrix transpositions (see article by DHB in Intl. J. of Supercomputer
!   Applications, Spring 1988, p. 82 - 87). This is the radix 2 version.
!   X is both the input and the output array, while Y is a scratch array.
!   Both X and Y must be dimensioned with 2 * N real cells, where N = 2^M.
!   The data in X are assumed to have real and imaginary parts separated
!   by N cells.  Before calling CFFTZ to perform an FFT, the array U must be
!   initialized by calling CFFTZ with IS set to 0 and M set to MX, where MX is
!   the maximum value of M for any subsequent call.  U must be dimensioned
!   with at least 2 * NX real cells, where NX = 2^MX.  M must be at least two.

!   David H. Bailey     May 4, 1988
*/

{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    int n = (int) pow(2, m);
    if (0 == is) {
	
	//!   Initialize the U array with sines and cosines in a manner that permits
	//!   stride one access at each FFT iteration.

	int nu = n;
	u[1-1] = m;
	int ku = 2;
	int kn = ku+nu;
	int ln = 1;

	for (int j = 1; j <= m; j++) {
	    mp_real t = mppic/ln;

	    //!   This loop is vectorizable.
	    for (int i = 0; i < ln; i++) {
		mp_real ti = i*t;
		mpcssnf(ti, u[i+ku-1], u[i+kn-1]);
	    }

	    ku += ln;
	    kn = ku+nu;
	    ln *= 2;
	}
	return;
    }

    //!   Check if input parameters are invalid.
    int mx = INT(u[1-1]);
    if ((is != 1 && is != -1) || m < 2 || m > mx) {
	cout << "CFFTZ: Either U has not been initialized, or else\n";
	cout << "one of the input parameters is invalid ";
	cout << is << m << mx << '\n';
	exit(1);
    }

/*
!   A normal call to CFFTZ starts here.  M1 is the number of the first variant
!   radix-2 Stockham iterations to be performed.  The second variant is faster
!   on most computers after the first few iterations, since in the second
!   variant it is not necessary to access roots of unity in the inner DO loop.
!   Thus it is most efficient to limit M1 to some value.  For many vector
!   computers, the optimum limit of M1 is 6.  For scalar systems, M1 should
!   probably be limited to 2.
*/
    int m1 = min(m/2, 2);
    // m1 = min(m/2, 6);
    int m2 = m-m1;
    int n2 = (int) pow(2, m1);
    int n1 = (int) pow(2, m2);
    int i;

    //!   Perform one variant of the Stockham FFT.
    for (int l = 1; l <= m1; l += 2) {
	fftz1(is, l, m, u, x, y);
	if (l == m1) goto lab140;
	fftz1(is, l+1, m, u, y, x);
    }

    //!   Perform a transposition of X treated as a N2 x N1 x 2 matrix.
    trans(n1, n2, x, y);

    //!   Perform second variant of the Stockham FFT from Y to X and X to Y.
    for (l = m1+1; l <= m; l += 2) {
	fftz2(is, l, m, u, y, x);
	if (l == m) goto lab180;
	fftz2(is, l+1, m, u, x, y);
    }

    goto lab160;
    //!   Perform a transposition of Y treated as a N2 x N1 x 2 matrix.
  lab140:;
    trans(n1, n2, y, x);

    //!   Perform second variant of the Stockham FFT from X to Y and Y to X.
    for (l = m1+1; l <= m; l += 2) {
	fftz2(is, l, m, u, x, y);
	if (l == m) goto lab160;
	fftz2(is, l+1, m, u, y, x);
    }
    goto lab180;

    //!   Copy Y to X.
  lab160:
    for (i = 0; i < 2*n; i++)
      x[i] = y[i];

  lab180:;
}

void fftz1(const int is, const int l, const int m, const mp_real* u,
	   mp_real* x, mp_real* y)

//!   Performs the L-th iteration of the first variant of the Stockham FFT.

{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    //!   Set initial parameters.

    int n = (int) pow(2, m);
    int mx = INT(u[1-1]);
    int nu = (int) pow(2, mx);
    int n1 = n/2;
    int lk = (int) pow(2, l-1);
    int li = (int) pow(2, m-l);
    int lj = 2*li;
    int ku = li+1;
    int kn = ku+nu;

    for (int k = 0; k < lk; k++) {
	int i11 = k*lj+1;
	int i12 = i11+li;
	int i21 = k*li+1;
	int i22 = i21+n1;

	//!   This loop is vectorizable.
	for (int i = 0; i < li; i++) {
	    mp_real u1 = u[ku+i-1];
	    mp_real u2 = is*u[kn+i-1];
	    mp_real x11 = x[i11+i-1];
	    mp_real x12 = x[i11+i+n-1];
	    mp_real x21 = x[i12+i-1];
	    mp_real x22 = x[i12+i+n-1];
	    mp_real t1 = x11-x21;
	    mp_real t2 = x12-x22;
	    y[i21+i-1] = x11+x21;
	    y[i21+i+n-1] = x12+x22;
	    y[i22+i-1] = u1*t1-u2*t2;
	    y[i22+i+n-1] = u1*t2+u2*t1;
	}
    }
}

void fftz2(const int is, const int l, const int m, const mp_real* u,
	   mp_real* x, mp_real* y)

//!   Performs the L-th iteration of the second variant of the Stockham FFT.
  
{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    //!   Set initial parameters.

    int n = (int) pow(2, m);
    int mx = INT(u[1-1]);
    int nu = (int) pow(2, mx);
    int n1 = n/2;
    int lk = (int) pow(2, l-1);
    int li = (int) pow(2, m-l);
    int lj = 2*lk;
    int ku = li+1;

    for (int i = 0; i < li; i++) {
	int i11 = i*lk+1;
	int i12 = i11+n1;
	int i21 = i*lj+1;
	int i22 = i21+lk;
	mp_real u1 = u[ku+i-1];
	mp_real u2 = is*u[ku+i+nu-1];

	//!   This loop is vectorizable.
	for (int k = 0; k < lk; k++) {
	    mp_real x11 = x[i11+k-1];
	    mp_real x12 = x[i11+k+n-1];
	    mp_real x21 = x[i12+k-1];
	    mp_real x22 = x[i12+k+n-1];
	    mp_real t1 = x11-x21;
	    mp_real t2 = x12-x22;
	    y[i21+k-1] = x11+x21;
	    y[i21+k+n-1] = x12+x22;
	    y[i22+k-1] = u1*t1-u2*t2;
	    y[i22+k+n-1] = u1*t2+u2*t1;
	}
    }
}

void trans(const int n1, const int n2, const mp_real* x, mp_real* y)

/*
!   Performs a transpose of the vector X, returning the result in Y.  X is
!   treated as a N1 x N2 complex matrix, and Y is treated as a N2 x N1 complex
!   matrix.  The complex data is assumed stored with real and imaginary parts
!   separated by N1 x N2 locations.
*/

{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)
    const int na = 32;
    const int nc = 32;
    int i, j;
    //DIMENSION X(2*N1*N2), Y(2*N1*N2), Z(NC,2*NC)

    int n = n1*n2;
    if (n1 >= n2)
      goto lab110;
    else
      goto lab130;

    //!   Scheme 1:  Perform a simple transpose in the usual way.
  lab110:
    for (j = 0; j < n2; j++) {
	int j1 = j+1;
	int j2 = j*n1+1;

	//!   This loop is vectorizable.

	for (i = 0; i < n1; i++) {
	    y[i*n2+j1-1] = x[i+j2-1];
	    y[i*n2+j1+n-1] = x[i+j2+n-1];
	}
    }
    goto lab260;

    //!   Scheme 2:  Perform a simple transpose with the loops reversed.
  lab130:
    for (i = 0; i < n1; i++) {
	int i1 = i*n2+1;
	int i2 = i+1;
	
	//!   This loop is vectorizable.

	for (j = 0; j < n2; j++) {
	    y[j+i1-1] = x[j*n1+i2-1];
	    y[j+i1+n-1] = x[j*n1+i2+n-1];
	}
    }
  lab260: return;
}
