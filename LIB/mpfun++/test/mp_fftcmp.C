/*
!   Performs real-to-complex and complex-to-real FFT routines using routines
!   RCFFTC, CRFFTC and CFFTC.  These routines use double complex data.

!   David H. Bailey     June 6, 1994

!   Translated into C++.
!   Siddhartha Chatterjee (RIACS)
!   18 August 1994
*/

#define DEBUG 1

#include <iostream.h>
#include "mpmod.H"
extern "C" {
    #include <time.h>
    #include <stdlib.h>
    #include <math.h>
}

double randlc(double&, const double);
void vranl(const int, double&, const double, double*);
void rcfftc(const int, const int, mp_complex*, mp_real*,
	    mp_complex*, mp_complex*);
void crfftc(const int, const int, mp_complex*, mp_complex*,
	    mp_complex*, mp_real*);
void cfftc(const int, const int, mp_complex*, mp_complex*, mp_complex*);
void fftc1(const int, const int, const int, const mp_complex*,
	   mp_complex*, mp_complex*);
void fftc2(const int, const int, const int, const mp_complex*,
	   mp_complex*, mp_complex*);
void trans(const int, const int, const mp_complex*, mp_complex*);

int min(const int m, const int n) {return m < n ? m : n;}

int main()

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    const int 		m = 8;
//    const int 		n = (int) pow(2, m);
    const int 		n = 256;
    const int 		it = 10;
    const double	aa = 1220703125.0;
    const double 	ss = 314159265.0;
    mp_real 		x[n+2], rn, se;
    mp_complex 		u[n], y[n/2+1], z[n/2+1];
    double 		s[n];
    time_t		t0, t1;

    //!   Initialize.  The SECOND function is assumed to be the timing function
    //!   on the given computer.  Change here and below if another name is used.

    mp_init();
    rn = mp_real(1.0)/n;
    double rr = 0.0;
    double xx;
    vranl(0, xx, aa, s);
    xx = ss;
    cfftc(0, m, u, y, z);
    time(&t0);

    //!  Perform IT iterations.
    for (int ii = 1; ii <= it; ii++) {
	vranl(n, xx, aa, s);
	for (int i = 0; i < n; i++) x[i] = s[i];

	rcfftc(-1, m, u, x, y, z);
	crfftc(1, m, u, z, y, x);
	se = 0.0;

	for (i = 0; i < n; i++) {
	    se += power(rn*x[i]-s[i], 2);
	}

	se = sqrt(se/n);
	cout << "Iteration " << ii << "   ";
	cout << "Error =\n" << se;
	if (se > mpeps) {
	    cout << "*** Error detected:\n";
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

double aint(double x)

{
    if (fabs(x) < 1.0) return 0.0;
    return floor(fabs(x))*(x >= 0 ? 1.0 : -1.0);
}

#if	0
double randlc(double& x, const double a)

/*
!   This routine returns a uniform pseudorandom double precision number in the
!   range (0, 1) by using the linear congruential generator

!   x_{k+1} = a x_k  (mod 2^46)

!   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
!   before repeating.  The argument A is the same as 'a' in the above formula,
!   and X is the same as x_0.  A and X must be odd double precision integers
!   in the range (1, 2^46).  The returned value RANDLC is normalized to be
!   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
!   the new seed x_1, so that subsequent calls to RANDLC using the same
!   arguments will generate a continuous sequence.

!   This routine should produce the same results on any computer with at least
!   48 mantissa bits in double precision floating point data.  On Cray systems,
!   double precision should be disabled.

!   David H. Bailey     October 26, 1990
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    static int ks = 0;
    static double r23, r46, t23, t46;

/*
!   If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
!   T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
!   by merely using the ** operator, in order to insure that the results are
!   exact on all systems.  This code assumes that 0.5D0 is represented exactly.
*/
    if (ks == 0) {
	r23 = r46 = t23 = t46 = 1.0;
	for (int i = 0; i < 23; i++) {
	    r23 *= 0.5;
	    t23 *= 2.0;
	}
	for (i = 0; i < 46; i++) {
	    r46 *= 0.5;
	    t46 *= 2.0;
	}
	ks = 1;
    }

    //!   Break A into two parts such that A = 2^23 * A1 + A2 and set
    //  X = N.
    double t1 = r23*a;
    double a1 = aint(t1);
    double a2 = a-t23*a1;

    //!   Break X into two parts such that X = 2^23 * X1 + X2, compute
    //!   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    //!   X = 2^23 * Z + A2 * X2  (mod 2^46).

    t1 = r23*x;
    double x1 = aint(t1);
    double x2 = x-t23*x1;
    t1 = a1*x2+a2*x1;
    double t2 = aint(r23*t1);
    double z = t1-t23*t2;
    double t3 = t23*z+a2*x2;
    double t4 = aint(r46*t3);
    x = t3-t46*t4;
    return r46*x;
}
#endif

void vranl(const int n, double& x, const double a, double *y)

/*
!   This routine generates N uniform pseudorandom double precision numbers in
!   the range (0, 1) by using the linear congruential generator

!   x_{k+1} = a x_k  (mod 2^46)

!   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
!   before repeating.  The argument A is the same as 'a' in the above formula,
!   and X is the same as x_0.  A and X must be odd double precision integers
!   in the range (1, 2^46).  The N results are placed in Y and are normalized
!   to be between 0 and 1.  X is updated to contain the new seed, so that
!   subsequent calls to VRANL using the same arguments will generate a
!   continuous sequence.  If N is zero, only initialization is performed, and
!   the variables X, A and Y are ignored.

!   This routine is the standard version designed for scalar or RISC systems.
!   However, it should produce the same results on any single processor
!   computer with at least 48 mantissa bits in double precision floating point
!   data.  On Cray systems, double precision should be disabled.

!   David H. Bailey     October 26, 1990
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    static int ks = 0;
    static double r23, r46, t23, t46;

    //!   If this is the first call to VRANL, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
    //!   T23 = 2 ^ 23, and T46 = 2 ^ 46.  See comments in RANDLC.
    if (0 == ks || 0 == n) {
	r23 = r46 = t23 = t46 = 1.0;
	for (int i = 0; i < 23; i++) {
	    r23 *= 0.5;
	    t23 *= 2.0;
	}
	for (i = 0; i < 46; i++) {
	    r46 *= 0.5;
	    t46 *= 2.0;
	}

	ks = 1;
	if (n == 0) return;
    }

    //!   Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.

    double t1 = r23*a;
    double a1 = aint(t1);
    double a2 = a-t23*a1;

    //!   Generate N results.   This loop is not vectorizable.
    for (int i = 0; i < n; i++) {
	//!   Break X into two parts such that X = 2^23 * X1 + X2, compute
	//!   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
	//!   X = 2^23 * Z + A2 * X2  (mod 2^46).

	t1 = r23*x;
	double x1 = aint(t1);
	double x2 = x-t23*x1;
	t1 = a1*x2+a2*x1;
	double t2 = aint(r23*t1);
	double z = t1-t23*t2;
	double t3 = t23*z+a2*x2;
	double t4 = aint(r46*t3);
	x = t3-t46*t4;
	y[i] = r46*x;
    }
}

void rcfftc(const int is, const int m, mp_complex *u, mp_real *x,
	    mp_complex *y, mp_complex *z)

/*
!   Computes the N-point real-to-complex FFT of X, where N = 2^M.   This 
!   uses an algorithm due to Swarztrauber, coupled with some fast methods of
!   methods of performing power-of-two matrix transforms.  X is the input
!   array, Y is a scratch array, and Z is the output array.  X must be
!   dimensioned with N + 2 real cells.  Y, and Z must be dimensioned with N/2+1
!   complex cells.  If desired, X and Y or X and Z may be the same array in
!   the calling program.  Before calling RCFFTC to perform an FFT, the array U
!   must be initialized by calling CFFTC with IS set to 0 and M set to MX,
!   where MX is the maximum value of M for any subsequent call.  U must be
!   dimensioned with at least 2 * N real cells.  M must be at least two.
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    mp_complex c0, c1, y0, y1;

    //!   Set initial parameters.
    int n = (int) pow(2, m);
    int mx = INT(u[0]);
    int n2 = n/2;
    int n21 = n2+1;
    int n22 = n2+2;
    int n4 = n/4;
    int n41 = n4+1;
    int ku = n/2+1;

    //!   Check if input parameters are invalid.
    if ((is != 1 && is != -1) || m < 3 || m > mx) {
	cout << "RCFFTC: Either U has not been initialized, or else ";
	cout << "one of the input parameters is invalid ";
	cout << is << n << mx << '\n';
	exit(1);
    }

    //!   Copy X to Y such that Y(k) = X(2k-1) + i X(2k).
    for (int k = 1; k <= n2; k++)
      y[k-1] = mp_complex(x[2*k-1-1], x[2*k-1]);

    //!   Perform a normal N/2-point FFT on Y.
    cfftc(is, m-1, u, y, z);

    //!   Reconstruct the FFT of X.
    z[1-1] = 2.0*(MP_REAL(y[1-1]) + aimag(y[1-1]));
    z[n41-1] = 2.0*mp_complex(MP_REAL(y[n41-1]), is*aimag(y[n41-1]));
    z[n21-1] = 2.0*(MP_REAL(y[1-1]) - aimag(y[1-1]));

    if (1 == is) {
	for (k = 2; k <= n4; k++) {
	    y0 = y[k-1];
	    y1 = conjg(y[n22-k-1]);
	    c0 = y0+y1;
	    c1 = u[ku+k-1]*(y0-y1);
	    z[k-1] = c0+c1;
	    z[n22-k-1] = conjg(c0-c1);
	}
    }
    else {
	for (k = 2; k <= n4; k++) {
	    y0 = y[k-1];
	    y1 = conjg(y[n22-k-1]);
	    c0 = y0+y1;
	    c1 = conjg(u[ku+k-1])*(y0-y1);
	    z[k-1] = c0+c1;
	    z[n22-k-1] = conjg(c0-c1);
	}
    }
}

void crfftc(const int is, const int m, mp_complex *u, mp_complex *x,
	    mp_complex *y, mp_real *z)

/*
!   Computes the N-point complex-to-real FFT of X, where N = 2^M.   This 
!   uses an algorithm due to Swarztrauber, coupled with some fast methods of
!   methods of performing power-of-two matrix transforms.  X is the input
!   array, Y is a scratch array, and Z is the output array.  X and Y must be
!   dimensioned with N/2+1 complex cells.  Z must be dimensioned with N + 2
!   real cells.  If desired, X and Y or X and Z may be the same array in
!   the calling program.  Before calling CRFFTC to perform an FFT, the array U
!   must be initialized by calling CFFTC with IS set to 0 and M set to MX,
!   where MX is the maximum value of M for any subsequent call.  U must be
!   dimensioned with at least 2 * N real cells.  M must be at least two.
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    mp_complex c0, c1, x0, x1;

    //!   Set initial parameters.
    int n = (int) pow(2,m);
    int mx = INT(u[0]);
    int n2 = n/2;
    int n21 = n2+1;
    int n22 = n2+2;
    int n4 = n/4;
    int n41 = n4+1;
    int ku = n/2+1;

    //!   Check if input parameters are invalid.
    if ((is != 1 && is != -1) || m < 3 || m > mx) {
	cout << "CRFFTC: Either U has not been initialized, or else ";
	cout << "one of the input parameters is invalid ";
	cout << is << n << mx << '\n';
	exit(1);
    }

    //!   Construct the input to CFFTC.
    y[1-1] = 0.5*mp_complex(MP_REAL(x[1-1])+MP_REAL(x[n21-1]),
      MP_REAL(x[1-1]) - MP_REAL(x[n21-1]));
    y[n4+1-1] = mp_complex(MP_REAL(x[n41-1]), -is*aimag(x[n41-1]));

    if (1 == is) {
	for (int k = 2; k <= n4; k++) {
	    x0 = x[k-1];
	    x1 = conjg(x[n22-k-1]);
	    c0 = x0+x1;
	    c1 = u[ku+k-1]*(x0-x1);
	    y[k-1] = 0.5*(c0+c1);
	    y[n22-k-1] = 0.5*conjg(c0-c1);
	}
    }
    else {
	for (int k = 2; k <= n4; k++) {
	    x0 = x[k-1];
	    x1 = conjg(x[n22-k-1]);
	    c0 = x0+x1;
	    c1 = conjg(u[ku+k-1])*(x0-x1);
	    y[k-1] = 0.5*(c0+c1);
	    y[n22-k-1] = 0.5*conjg(c0-c1);
	}
    }

    //!   Perform a normal N/2-point FFT on Y.
    cfftc(is, m-1, u, y, x);

    //!   Copy Y to Z such that Y(k) = Z(2k-1) + i Z(2k).
    for (int k = 1; k <= n2; k++) {
	z[2*k-1-1] = MP_REAL(y[k-1]);
	z[2*k-1] = aimag(y[k-1]);
    }
}

void cfftc(const int is, const int m, mp_complex *u, mp_complex *x,
	   mp_complex *y)

/*
!   Computes the N-point complex-to-complex FFT of X, where N = 2^M.   This 
!   uses an algorithm due to Swarztrauber, coupled with some fast methods of
!   methods of performing power-of-two matrix transforms.  X is the input
!   array, Y is a scratch array, and Z is the output array.  X, Y, and Z must
!   be dimensioned with N complex cells.  If desired, X and Y or X and Z may
!   be the same array in the calling program.  Before calling CFFTC to perform
!   an FFT, the array U must be initialized by calling CFFTC with IS set to 0
!   and M set to MX, where MX is the maximum value of M for any subsequent
!   call.  U must be dimensioned with at least 2 * N real cells.  M must be at
!   least two.
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    mp_real t, ti, t1, t2;
//    int n = (int) pow(2, m);
    const int n = 256;
    mp_complex z[n];
    int i;

    if (0 == is) {

	//!   Initialize the U array with sines and cosines in a manner that permits
	//!   stride one access at each radix-4 FFT iteration.  If M is odd, then the
	//!   final set of results are for a radix-2 FFT iteration.

	u[1-1] = m;
	int ku = 2;
	int ln = 1;

	for (int j = 1; j <= m; j++) {
	    t = mppic/ln;
	    for (i = 0; i < ln; i++) {
		ti = i*t;
		mpcssnf(ti, t1, t2);
		u[i+ku-1] = mp_complex(t1, t2);
	    }

	    ku += ln;
	    ln *= 2;
	}
	return;
    }

    //!   Check if input parameters are invalid.
    int mx = INT(u[1-1]);
    if ((is != 1 && is != -1) || m < 2 || m > mx) {
	cout << "CFFTC: Either U has not been initialized, or else ";
	cout << "one of th input parameters is invalid ";
	cout << is << n << mx << '\n';
	abort();
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
    int m2 = m-m1;
    int n2 = (int) pow(2, m1);
    int n1 = (int) pow(2, m2);

    //!   Perform the first variant of the radix-2 Stockham FFT.
    for (int l = 1; l <= m1; l += 2) {
	fftc1(is, l, m, u, x, y);
	if (l == m1) goto lab150;
	fftc1(is, l+1, m, u, y, x);
    }

    goto lab170;
  lab150:;
    trans(n1, n2, y, x);

    //!   Perform the second variant of the radix-2 Stockham FFT.
    for (l = m1+1; l <= m; l += 2) {
	fftc2(is, l, m, u, x, y);
	if (l == m) goto lab190;
	fftc2(is, l+1, m, u, y, x);
    }

    goto lab210;
  lab170:;
    trans(n1, n2, x, y);

    //!   Perform the second variant of the radix-2 Stockham FFT.
    for (l = m1+1; l <= m; l += 2) {
	fftc2(is, l, m, u, y, x);
	if (l == m) goto lab210;
	fftc2(is, l+1, m, u, x, y);
    }

    //!   Copy Y to X if the last operation above stored into Y.
  lab190:
    for (i = 0; i < n; i++) x[i] = y[i];
  lab210:
    return;
}

void fftc1(const int is, const int l, const int m, const mp_complex
	   *u, mp_complex *x, mp_complex *y)

/*
!   Performs the L-th iteration of the first variant of the radix-2 Stockham
!   FFT.  Complex version.
*/

{
    mp_complex uj, x1, x2;

    //!   Set initial parameters.
    int n = (int) pow(2, m);
    int ni = n/2;
    int lk = (int) pow(2, l-1);
    int lj = ni/lk;
    int ku = lj+1;

    for (int k = 0; k < lk; k++) {
	int k11 = 2*k*lj+1;
	int k12 = k11+lj;
	int k21 = k*lj+1;
	int k22 = k21+ni;

	if (1 == is) {
	    for (int j = 0; j < lj; j++) {
		uj = u[ku+j-1];
		x1 = x[k11+j-1];
		x2 = x[k12+j-1];
		y[k21+j-1] = x1+x2;
		y[k22+j-1] = uj*(x1-x2);
	    }
	}
	else {
	    for (int j = 0; j < lj; j++) {
		uj = conjg(u[ku+j-1]);
		x1 = x[k11+j-1];
		x2 = x[k12+j-1];
		y[k21+j-1] = x1+x2;
		y[k22+j-1] = uj*(x1-x2);
	    }
	}
    }
}

void fftc2(const int is, const int l, const int m, const mp_complex
	   *u, mp_complex *x, mp_complex *y)

/*
!   Performs the L-th iteration of the second variant of the radix-2 Stockham
!   FFT.  Complex version.
*/

{
    mp_complex uj, x1, x2;

    //!   Set initial parameters.
    int n = (int) pow(2, m);
    int ni = n/2;
    int lk = (int) pow(2, l-1);
    int lj = ni/lk;
    int ku = lj+1;

    for (int j = 0; j < lj; j++) {
	int j11 = j*lk+1;
	int j12 = j11+ni;
	int j21 = 2*j*lk+1;
	int j22 = j21+lk;
	uj = u[ku+j-1];
	if (-1 == is) uj = conjg(uj);

	for (int k = 0; k < lk; k++) {
	    x1 = x[j11+k-1];
	    x2 = x[j12+k-1];
	    y[j21+k-1] = x1+x2;
	    y[j22+k-1] = uj*(x1-x2);
	}
    }
}

void trans(const int n1, const int n2, const mp_complex *x, mp_complex *y)

/*
!   Performs a transpose of the vector X, returning the result in Y.  X is
!   treated as a N1 x N2 complex matrix, and Y is treated as a N2 x N1 complex
!   matrix.  The complex data is assumed stored with real and imaginary parts
!   interleaved as in the COMPLEX type.

!   David H. Bailey      April 28, 1987
*/

{
    int n = n1*n2;
    int i, j;
/*
!   Perform one of three techniques, depending on N.  The best strategy varies
!   with the computer system.  This strategy is best for scalar systems.
*/

    if (n1 >= n2)
      goto lab100;
    else
      goto lab120;

  lab100:
    for (j = 0; j < n2; j++) {
	for (i = 0; i < n1; i++) {
	    y[i*n2+j+1-1] = x[i+j*n1+1-1];
	}
    }
    goto lab200;

/*
!   Scheme 2:  Perform a simple transform with the loops reversed.  This is
!   usually the best on vector computers if N1 is odd, or if both N1 and N2 are
!   small, and N2 is larger than N1.
*/

  lab120:
    for (i = 0; i < n1; i++) {
	for (j = 0; j < n2; j++) {
	    y[j+i*n2+1-1] = x[j*n1+i+1-1];
	}
    }
  lab200:;
}
