#include <iostream.h>
#include <fstream.h>
#include "mpmod.H"
extern "C" {
    #include <assert.h>
    #include <time.h>
    #include <stdlib.h>
    #include <stdio.h>
    #include <math.h>
}

extern void mppslq(const int, const double, const int, const int,
		   mp_real*, int&, mp_real*);
extern void initmp(const int, const int, mp_real**, mp_real**,
		   mp_real**, mp_real*, mp_real*, mp_real*);
extern void itermp(const int, const double, const int, const int,
		   mp_real**, mp_real**, mp_real**, mp_real*, int&);
extern double bound(const int, mp_real**);
extern void matout(const int, const int, mp_real**);
extern void matout(const int, const int, const int, mp_real**);
extern void matout(const int, const int, mp_real*);
extern double drand();

int min(const int x, const int y) {return x < y ? x : y;}

double min(const double x, const double y) {return x < y ? x : y;}

int max(const int x, const int y) {return x > y ? x : y;}

double max(const double x, const double y) {return x > y ? x : y;}

double aint(double x)

{
    if (fabs(x) < 1.0) return 0.0;
    return floor(fabs(x))*(x >= 0 ? 1.0 : -1.0);
}

double sign(const double a, const double b)
{return fabs(a)*(b >= 0 ? 1 : -1);}

double anint(double a) {return a > 0 ? aint(a+0.5) : aint(a-0.5);}

int main()

/*
!  Test program for MPPSLQ, a new version of Helaman Ferguson's PSLQ algorithm
!  for detecting integer relations in a vector of real numbers.
!  One-level Fortran-90 MP version.

!  David H. Bailey     June 3, 1994

!  Converted into C++.
!  Siddhartha Chatterjee (RIACS)
!  27 August 1994.

!  These parameters are set in the PARAMETER statement below.  Some other
!  parameters are set in other routines.
!    IDB = Debug level (0 -- 3).
!    ITR = Number of random trials when KQ = 1 or 2.
!    GAM = Gamma, the PSLQ algorithm parameter.  Typically set to either
!          Sqrt (4/3) = 1.154700538 or 4/3 = 1.333333333.
!    N   = Integer relation vector length.
!    KQ  = 0 for the algebraic case [1, AL, AL^2, ... AL^(N-1)], where
!          AL = 3^(1/KR) - 2^(1/KS).  Set N = KR * KS + 1 and ITR = 1.
!        = 1 for a random integer vector with a relation.
!        = 2 for a random integer vector without a relation.
!        = 3 for testing algebraic relations of a number read from a file.
!        = 4 for testing additive relations of numbers read from a file.
!        = 5 for testing multiplicative relations of numbers read from a file.
!        = 6 for custom input.
!    KR  = Degree of root of 3 when KQ = 0.
!    KS  = Degree of root of 2 when KQ = 0.
!    IDV = Number of digits in the random vector values (except the last)
!          when KQ = 1 or 2.
!    MR  = Log_10 of the max. random relation coefficients when KQ = 1 or 2.
!    MN  = Log_10 of the max. acceptable norm of relations.  If KQ = 1, set
!          MN = MR + 1.
!    RM  = Maximum size of random relation coefficients when KQ = 1 or 2
!          (may be set as 10^MR or to some other value).
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    const int idb = 2;
    const int itr = 1;
    const double gam = 1.154700538;
    const int n = 21, kq = 0;
    const int kr = 4;
    const int ks = 5;
    const int idv = 0;
    const int mr = 20;
    const int mn = mr+1;
    double rm = pow(10, mr);
    mp_real r[n], x[n], al, t1, ti;
    double dr[n], r1[n];

    mp_init();
    int ndp = mpipl;
    cout << "PSLQ Test Program\n";
    cout << "No. trials = " << itr << "   ";
    cout << "Gamma = " << gam << "   ";
    cout << "N = " << n << "   ";
    cout << "KQ = " << kq << "   ";
    cout << "KR = " << kr << "   ";
    cout << "KS = " << ks << '\n';
    cout << "   IDV = " << idv << "   ";
    cout << "MR = " << mr << "   ";
    cout << "MN = " << mn << "   ";
    cout << "RM = " << rm << '\n';
    cout << "Full precision = " << ndp << " digits\n";

    double tm = 0.0;
    int kps = 0;
    int kcs = 0;
    double d1;
    int i;
//    FILE* fp;

    for (int k = 1; k <= itr; k++) {
	cout << "Start trial " << k << '\n';
	switch (kq) {
	  case 0:
	    //!  This code generates AL = 3^(1/KR) - 2^(1/KS).  AL is algebraic of degree
	    //!  KR * KS.  Use N = KR * KS + 1 to recover the polynomial satisfied by AL.

	    al = mpnrtf(mp_real(3.0), kr) - mpnrtf(mp_real(2.0), ks);
	    break;
	  case 1:
	    //!  Construct a random input vector X with a relation.
	    d1 = 2.0*rm+1.0;
	    ti = power(mp_real(10.0), idv);
	    t1 = 0.0;
	    
	    for (i = 1; i < n; i++) {
		r1[i-1] = aint(d1*drand())-rm;
		x[i-1] = aint(t1*mp_rand());
		t1 += r1[i-1]*x[i-1];
	    }

	    r1[n-1] = dble(-sign(mp_real(1.0), t1));
	    x[n-1] = abs(t1);
	    cout << "Constructed relation:\n";
	    for (i = 0; i < n; i++)
	      cout << r1[i] << '\n';
	    break;
	  case 2:
	    ti = power(mp_real(10.0), idv);

	    //!  Construct a random input vector X without a relation.

	    for (i = 0; i < n; i++)
	      x[i] = aint(t1*mp_rand());
	    break;
	  case 3:
	    //!  Read an algebraic constant from a file.

//	    fp = fopen("pslq.inp", "r");
//	    fp >> al;
	    {
		ifstream from("pslq.inp");
		from >> al;
	    }
	    break;
	  case 4:
	    //!  Read constants from a file for additive test.

/*	    fp = fopen("pslq.inp", "r");
	    for (i = 0; i < n; i++) {
		fp >> al;
		x[i] = al;
	    }*/
	    {
		ifstream from("pslq.inp");
		for (i = 0; i < n; i++) {
		    from >> al;
		    x[i] = al;
		}
	    }
	    break;
	  case 5:
	    //!  Read constants from a file for multiplicative test.

/*	    fp = fopen("pslq.inp", "r");
	    for (i = 0; i < n; i++) {
		fp >> al;
		x[i] = log(al);
	    }*/
	    {
		ifstream from("pslq.inp");
		for (i = 0; i < n; i++) {
		    from >> al;
		    x[i] = log(al);
		}
	    }
	    break;
	  case 6:
	    //!  Produce X vector by a custom scheme.
	    break;
	}

	//!  If KQ is 0 or 3, generate X = [1, AL, AL^2, ..., AL^(N-1)].

	if (0 == kq || 3 == kq) {
	    x[1-1] = 1.0;
	    for (int i = 2; i <= n; i++) {
		x[i-1] = al*x[i-1-1];
	    }
	}

	//!  Perform relation search.

	int iq = 0;
	time_t tm0, tm1;
	time(&tm0);
	mppslq(idb, gam, n, mn, x, iq, r);
	time(&tm1);

	//!  Check if original relation was recovered.

	if (iq != 0) {
	    tm += difftime(tm1, tm0);

	    for (int i = 1; i <= n; i++) {
		if (abs(r[i-1]) < 1e15)
		  dr[i-1] = dble(r[i-1]);
		else
		  dr[i-1] = dble(sign(mp_real(999999999999999.0), r[i-1]));
	    }

	    cout << "Recovered relation:\n";
	    for (i = 0; i < n; i++) cout << dr[i] << '\n';

	    if (kq >= 1) {
		double d1 = 0.0;
		double d2 = 0.0;

		for (i = 1; i <= n; i++) {
		    d1 += abs(dble(r[i-1])-r1[i-1]);
		    d2 += abs(dble(r[i-1])+r1[i-1]);
		}

		if (abs(d1) <= 1e-6 || abs(d2) <= 1e-6) {
		    kcs++;
		    cout << "Original relation recovered.\n";
		}
		else
		  kps++;
	    }
	}
//	cout << "Total CPU time = " << difftime(tm1, tm0) << '\n';
    }

    if (itr > 1) {
	cout << "No. partial successes = " << kps << '\n';
	cout << "No. complete successes = " << kcs << '\n';
	int kss = kps+kcs;
	if (kss != 0) {
	    tm /= kss;
	    cout << "Ave. CPU Time of PS or CS runs = " << tm << '\n';
	}
    }
    return 0;
}

void mppslq(const int idb, const double gam, const int n, const int
	    mn, mp_real *x, int& iq, mp_real *r)
/*
!  The following parameters are set in this routine:
!    IPI = Iteration print interval when IDB = 2.
!    ITM = Maximum iteration count.  Run is aborted when this is exceeded.
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    const int ipi = 100;
    const int itm = 100000;
    mp_real **a, **b, **h, **ss, *s, *w1;
    mp_real *w2, *y, ty1, ty2;
    int i, j;
    double dm;

    a = new mp_real*[n];
    b = new mp_real*[n];
    h = new mp_real*[n];
    ss = new mp_real*[n];
    s = new mp_real[n];
    w1 = new mp_real[n];
    w2 = new mp_real[n];
    y = new mp_real[n];
    for (int ii = 0; ii < n; ii++) {
	a[ii] = new mp_real[n];
	b[ii] = new mp_real[n];
	h[ii] = new mp_real[n-1];
	ss[ii] = new mp_real[n];
    }

    //!  Initialize.
    double tm2 = 0.0;
    iq = 0;
    int it = 0;
    double rn = 0.0;
    double rb = pow(10.0, mn);
    initmp(idb, n, a, b, h, s, x, y);

  lab100:
    it++;
    if (it > itm) {
	cout << "Iteration limit exceeded " << itm << '\n';
	goto lab170;
    }

    if (0 == it%ipi) {
	//!  Compute norm bound.

	double d1 = bound(n, h);
	rn = max(rn, d1);

	if (idb >= 2) {
	    ty1 = 1e100;
	    ty2 = 0.0;

	    for (i = 0; i < n; i++) {
		ty1 = min(ty1, abs(y[i]));
		ty2 = max(ty2, abs(y[i]));
	    }

	    cout << "Iteration " << it << '\n';
	    cout << "Norm bound = " << d1 << "    ";
	    cout << "Max. bound = " << rn << '\n';
	    cout << "Min, max of Y:\n";
	    cout << ty1 << ty2;
	}
	if (rn > rb) {
	    cout << "Norm bound limit exceeded " << rb << '\n';
	    goto lab170;
	}
    }

    int izm;

    //!  Perform MP reductions.
    itermp(idb, gam, it, n, a, b, h, y, izm);
    if (-2 == izm) goto lab170;
    if (izm > 0) goto lab130;
    goto lab100;

    //!  A relation has been detected.  Output the final norm bound.
  lab130:
    if (idb >= 1) {
	cout << "Relation detected.\n";
	cout << "No. iterations = " << it << '\n';
	cout << "Max. bound = " << rn << '\n';
    }
    dm = 1e100;

    //!  If there is more than one relation, select the one with smallest norm.

    for (j = 1; j <= n; j++) {
	if (abs(y[j-1]) <= mpeps) {
	    double d1 = 0.0;
	    for (i = 1; i <= n; i++)
	      d1 += pow(dble(b[i-1][j-1]), 2);
	    d1 = sqrt(d1);
	    if (idb >= 1) {
		cout << "Index of relation = " << j << "    ";
		cout << "Norm = " << d1 << '\n';
	    }
	    if (idb >= 2) {
		cout << "Residual = \n" << y[j-1];
		cout << "Relation:\n";
//		matout(1, n, b[1-1][j-1]);
		matout(1, n, j-1, b);
	    }
	    if (d1 < dm) {
		dm = d1;
		for (i = 1; i <= n; i++)
		  r[i-1] = b[i-1][j-1];
	    }
	}
    }

    if (dm <= rb)
      iq = 1;
    else {
	if (idb >= 2) cout << "Relation is too large.\n";
    }
  lab170:
    return;
}

void initmp(const int idb, const int n, mp_real **a, mp_real **b,
	    mp_real **h, mp_real *s, mp_real *x, mp_real *y)

//!  Initializes MP arrays at the beginning.
  
{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    if (idb >= 2) {
	cout << "Input X vector:\n";
	matout(1, n, x);
    }

    //!  Set A and B to the identity matrix.

    for (int j = 0; j < n; j++) {
	for (int i = 0; i < n; i++) {
	    a[i][j] = b[i][j] = 0.0;
	}
	a[j][j] = b[j][j] = 1.0;
    }

    mp_real t1 = 0.0;

    //!  Compute the S vector, the square root of the partial sum of squares of X,
    //!  and compute the Y vector, which is the normalized X vector.

    for (int i = n; i >= 1; i--) {
	t1 += power(x[i-1], 2);
	s[i-1] = sqrt(t1);
    }

    t1 = 1.0/s[1-1];

    for (i = 0; i < n; i++) {
	y[i] = t1*x[i];
	s[i] *= t1;
    }
    
    //!  Compute the initial H matrix.

    for (i = 1; i <= n; i++) {
	for (j = i+1; j <= n-1; j++)
	  h[i-1][j-1] = 0.0;

	if (i <= n-1)
	  h[i-1][i-1] = s[i+1-1]/s[i-1];

	for (j = 1; j <= i-1; j++)
	  h[i-1][j-1] = -y[i-1]*y[j-1]/(s[j-1]*s[j+1-1]);
    }

    //!  Perform full reduction on H, updating A and B.

    for (i = 2; i <= n; i++) {
	for (j = i-1; j >= 1; j--) {
	    if (abs(h[i-1][j-1]) > 0.5*abs(h[j-1][j-1])) {
		mp_real tm = anint(h[i-1][j-1]/h[j-1][j-1]);
		y[j-1] += tm*y[i-1];

		for (int k = 1; k <= j; k++)
		  h[i-1][k-1] -= tm*h[j-1][k-1];

		for (k = 1; k <= n; k++) {
		    a[i-1][k-1] -= tm*a[j-1][k-1];
		    b[k-1][j-1] += tm*b[k-1][i-1];
		}
	    }
	}
    }

    if (idb >= 3) {
	cout << "Initial H:\n";
	matout(n, n-1, h);
	cout << "Initial A matrix\n";
	matout(n, n, a);
	cout << "Initial B matrix\n";
	matout(n, n, b);
	cout << "Initial Y:\n";
	matout(1, n, y);
    }
}

void itermp(const int idb, const double gam, const int it, const int
	    n, mp_real **a, mp_real **b, mp_real **h, mp_real *y, int&
	    izm)

//!  This performs one iteration of the PSLQ algorithm using MP arithmetic.

{
    //IMPLICIT TYPE (MP_REAL) (A-H, O-Z)

    if (idb >= 3)
      cout << "Iteration " << it << "   MP iteration\n";
    mp_real tl1 = 1.0/mpeps;

    //!  Select IM = I such that GAM^I * |H(I,I)| is maximal.

    izm = 0;
    mp_real t1 = 0.0, t2, t3, t4;
    int im;

    for (int i = 1; i < n; i++) {
	t2 = pow(gam, i) * abs(h[i-1][i-1]);
	if (t2 > t1) {
	    im = i;
	    t1 = t2;
	}
    }

    if (idb >= 3) cout << "IM = " << im << '\n';
    int im1 = im+1;

    //!  Perform block reduction on H, updating A and B.

    for (i = im1; i <= n; i++) {
	int jl = min(i-1, im1);

	for (int j = jl; j >= 1; j--) {
	    if (abs(h[i-1][j-1]) > 0.5*abs(h[j-1][j-1])) {
		mp_real tm = anint(h[i-1][j-1]/h[j-1][j-1]);
		y[j-1] += tm*y[i-1];

		for (int k = 1; k <= j; k++)
		  h[i-1][k-1] -= tm*h[j-1][k-1];

		for (k = 1; k <= n; k++) {
		    a[i-1][k-1] -= tm*a[j-1][k-1];
		    b[k-1][j-1] += tm*b[k-1][i-1];
		}
	    }
	}
    }

    //!  Find the max of the updated A and the min of Y.

    mp_real ta = 0.0;
    mp_real ty = 1e100;

    for (int j = 1; j <= im1; j++) {
	ty = min(ty, abs(y[j-1]));

	for (i = im1; i <= n; i++)
	  ta = max(ta, abs(a[i-1][j-1]));
    }

    if (ta > tl1) {
	if (idb >= 1) {
	    cout << "Iteration " << it << '\n';
	    cout << "ITERMP: Large entry in A:\n";
	    cout << ta;
	    cout << "Run aborted.\n";
	}
	izm = -2;
	return;
    }

    if (abs(h[n-1-1][n-1-1]) <= mpeps) {
	if (idb >= 2) {
	    cout << "Iteration " << it << '\n';
	    cout << "ITERMP: Small value in H:\n";
	    cout << h[n-1-1][n-1-1];
	}
	izm = 1;
    }

    if (ty <= mpeps) {
	if (idb >= 2) {
	    cout << "Iteration " << it << '\n';
	    cout << "ITERMP: Small value in Y:\n";
	    cout << ty;
	}
	izm = 1;
    }

    //!  Exchange the IM and IM+1 rows of A, the columns of B and the rows
    //!  of H.

    t1 = y[im-1];
    y[im-1] = y[im1-1];
    y[im1-1] = t1;

    for (i = 1; i <= n; i++) {
	t1 = a[im-1][i-1];
	a[im-1][i-1] = a[im1-1][i-1];
	a[im1-1][i-1] = t1;
	t1 = b[i-1][im-1];
	b[i-1][im-1] = b[i-1][im1-1];
	b[i-1][im1-1] = t1;
    }

    for (i = 1; i < n; i++) {
	t1 = h[im-1][i-1];
	h[im-1][i-1] = h[im1-1][i-1];
	h[im1-1][i-1] = t1;
    }

    //!  Update H with permutation produced above.

    if (im <= n-2) {
	t1 = h[im-1][im-1];
	t2 = h[im-1][im1-1];
	t3 = sqrt(power(t1,2)+power(t2,2));
	t1 /= t3;
	t2 /= t3;

	for (i = im; i <= n; i++) {
	    t3 = h[i-1][im-1];
	    t4 = h[i-1][im1-1];
	    h[i-1][im-1] = t1*t3+t2*t4;
	    h[i-1][im1-1] = -t2*t3+t1*t4;
	}
    }

    if (idb >= 3) {
	cout << "Updated A matrix:\n";
	matout(n, n, a);
	cout << "Updated B matrix:\n";
	matout(n, n, b);
	cout << "Updated H matrix:\n";
	matout(n, n-1, h);
	cout << "Updated Y:\n";
	matout(1, n, y);
    }
}

double bound(const int n, mp_real **h)

//!  This computes the norm bound.

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)

    double t1 = 0.0;

    for (int i = 1; i <= n; i++) {
	double t2 = 0.0;
	for (int j = 1; j < n; j++) {
	    t2 += pow(dble(h[i-1][j-1]), 2);
	}
	t1 = max(t1, t2);
    }
    return 1.0/sqrt(t1);
}

void matout(const int n1, const int n2, mp_real **a)

//!  This outputs the MP matrix A.

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)

    for (int i = 1; i <= n1; i++) {
	cout << "Row " << i << '\n';
	for (int j = 1; j <= n2; j++)
	  cout << a[i-1][j-1];
    }
}

void matout(const int n1, const int n2, const int c, mp_real **a)

//!  This outputs column c of the MP matrix A.

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)

    for (int i = 1; i <= n1; i++) {
	cout << "Row " << i << '\n';
	cout << a[i-1][c];
    }
}

void matout(const int n1, const int n2, mp_real *a)

//!  This outputs the MP matrix A.

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)

    assert(1 == n1 || 1 == n2);
    for (int i = 1; i <= n1; i++) {
	cout << "Row " << i << '\n';
	for (int j = 1; j <= n2; j++)
	  cout << (1 == n1 ? a[j-1] : a[i-1]);
    }
}

double drand()

/*
!  This routine returns a pseudorandom DP floating number nearly uniformly
!  distributed between 0 and 1 by means of a linear congruential scheme.
!  2^28 pseudorandom numbers with 30 bits each are returned before repeating.
*/

{
    //IMPLICIT DOUBLE PRECISION (A-H, O-Z)
    const double f7 = 78125.0;
    double r30 = pow(0.5, 30);
    double t30 = pow(2, 30);
    static double sd = 314159265.0;

    double t1 = f7*sd;
    double t2 = aint(r30*t1);
    sd = t1-t30*t2;
    return r30*sd;
}
