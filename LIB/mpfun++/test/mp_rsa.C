/*
!   Tests the performance and integrity of a supercomputer system by
!   encrypting and decrypting with the RSA algorithm.  Converted to 
!   Fortran-90 based multiprecision.

!   David H. Bailey     June 3, 1994
!
!   Translated into C++.
!
!   Siddhartha Chatterjee (RIACS)
!   9 August 1994
*/

#include "mpmod.H"
#include <iostream.h>
extern "C" {
    #include <time.h>
    #include <stdlib.h>
    #include <math.h>
}

int idivm(int idb, mp_integer& n, int *ipt, int kp)

{

    //!   Returns the first prime divisor of the multiprecision N.  Subroutine
    //!   PRIMES must have been previously called to initialize the array IPT.


    mp_integer k;
    int ret, i;

    //!   Check 3 separately (must be congruent to 2 mod 3).
    int id = 3;
    k = n/id;
    int l = INT(n - k*id);
    if (l == 0) goto lab120;
    if (l == 1) goto lab130;

    for (i = 2; i <= kp; i++) {
	id = ipt[i-1];
	k = n/id;
	l = INT(n-k*id);
	if (l == 0) goto lab120;
    }

    id = 0;

  lab120:
    ret = id;
    if (idb >= 3 && id != 0)
      cout << "N is not prime.  Divisible by " << id << '\n';
    goto lab200;

  lab130:
    ret = -3;
    if (idb >= 3) cout << "N is not congruent to 2 mod 3\n";

  lab200:
    return ret;
}

void ranprm(int idb, mp_integer& ip, int *ipt, int kp, int nd4)

{

    //!   Finds a large prime number.  Subroutine PRIMES must have
    //been previously 
    //!   called to initialize the array IP and determine KP.

    mp_integer id, iq, ix, iy, iz, k, n, n1;
    const int it = 10;
    int ii, jj, ik;

    int nt = 20*nd4;
    id = power(mp_integer(10), nd4)/2;

    for (int j = 1; j <= nt; j++) {

	//!   Generate a trial prime N that is odd.
//	n1 = 2*mp_integer(id*mp_rand());
	n1 = id*mp_rand();
	n1 = 2*n1;
	n = n1+1;
	if (idb >= 2) cout << "PRIME TRIAL " << j << '\n';
	if (idb >= 3) cout << "N = " << n;

	//!   Perform quick-and-dirty check for primality.
	int kd = idivm(idb, n, ipt, kp);
	if (kd != 0) goto lab200;

	//!   Quick test failed to find any divisors.  Proceed with
	//advanced test.
	iq = n1;

	for (ik = 0; ik <= 30; ik++) {
	    k = iq/2;
	    if (iq != 2*k) goto lab110;
	    iq = k;
	}

	//!   Now N = 1 + Q * 2^IK, where Q is odd.
      lab110:
	for (ii = 1; ii <= it; ii++) {
	    ix = n*mp_rand();
	    iy = 1;
	    iz = ix;
	  lab120:
	    k = iq/2;
	    if (iq != 2*k) iy = mod(iz*iy,n);
	    iq = k;
	    if (iq == 0) goto lab130;
	    iz = mod(iz*iz,n);
	    goto lab120;

	    //!   Now Y = X^Q (mod N).
	  lab130:
	    jj = 0;
	    if (iy == 1) goto lab150;

	    //!   Check primality conditions.
	  lab140:
	    jj++;
	    if (jj > ik) goto lab160;
	    if (iy == n1) goto lab150;
	    if (iy == 1) goto lab160;
	    iy = mod(iy*iy,n);
	    goto lab140;
	  lab150:;
	}

	//!   N is probably prime.
	if (idb >= 1) {
	    cout << "Prime found at trial " << j << '\n';
	    cout << n;
	}
	goto lab210;

	//!   N is definitely not prime.
      lab160:
	if (idb >= 2)
	  cout << "N is not prime.  No. advanced iterations = " << ii
	    << '\n';
      lab200:;
    }

    cout << "Failed to find a prime.  No. trials = " << nt << '\n';
    exit(1);

  lab210:
    ip = n;
}

int idiv(int n)

{
    //!   Returns the first prime divisor of N (single precision).

    static const int id[8] = {1, 7, 11, 13, 17, 19, 23, 29};
    int kd, i, k, l, m;

    if (n%3 == 0) {
	kd = 3;
	goto lab200;
    }

    if (n%5 == 0) {
	kd = 5;
	goto lab200;
    }

    kd = 0;
    if (n == 29) goto lab200;
    k = 2;
    l = 0;
    m = int(sqrt(double(n)));
  lab100:
    for (i = k; i <= 8; i++) {
	kd = id[i-1] + l;
	if (n%kd == 0) goto lab200;
    }

    kd = 0;
    l += 30;
    if (l > m) goto lab200;
    k = 1;
    goto lab100;

  lab200:
    return kd;
}

void primes(int idb, int *ipt, int kp)

{

    //!   Places the first KP prime numbers in IPT.
    static const int ipi[9] = {3, 5, 7, 11, 13, 17, 19, 23, 29};

    for (int i = 1-1; i < 9; i++)
      ipt[i] = ipi[i];

    int k = 9;

    for (i = 31; i <= 1000000; i += 2) {
	if (idiv(i) == 0) {
	    ipt[k] = i;
	    k++;
	    if (k == kp) return;
	}
    }

    cout << "Table IPT not filled: " << k<< '\n';
    exit(1);
}

void expm(int idb, mp_integer& ia, mp_integer& ib, mp_integer& ic,
	  mp_integer& ip)

{
    //!   Computes  IP = IA ** IB (mod IC)  for integers IA, IB, IC.

    mp_integer k, l, l2;

    k = ia;
    l = ib;
    ip = 1;
    int j = 0;


    //!   Perform binary algorithm for exponentiation modulo IC.
  lab100:
    j++;
    l2 = l/2;
    if (l != 2*l2) ip = mod(k*ip,ic);
    k = mod(k*k,ic);
    l = l2;
    if (l != 0) goto lab100;
    if (idb >= 2) cout << "EXPM iterations = " << j << '\n';
}

int main()

{
    mp_integer	ia, ib, id, ie, ip, iq, ix, iy, iz;
    time_t	tm0, tm1, tm2, tt0, tt1;
    const int	n1 = 10, n2 = 10, kp = 200, idb = 0;
    int		ipt[kp];

    //!  Initialize.
    mp_init();
    time(&tt0);
    int nd = mpipl-15;
    int nd4 = nd/4;
    primes(idb, ipt, kp);

    //!   Begin outer loop.
    for (int k = 1; k <= n2; k++) {
	time(&tm0);
	cout << "Outer trial " << k << '\n';

	//!   Step 1:  obtain two large primes, each congruent to 2
	//mod 3, with a size
	//!   somewhat less than half the machine precision level.
	ranprm(idb, ip, ipt, kp, nd4);
      lab100:;
	ranprm(idb, iq, ipt, kp, nd4);
	if (ip == iq) goto lab100;
	if (idb >= 1) {
	    cout << "p = " << ip;
	    cout << "q = " << iq;
	}

	//!   Step 2:  Compute  A = P * Q  and  B = (P-1) * (Q-1).
	ia = ip*iq;
	ib = (ip-1)*(iq-1);

	//!   Step 3:  Set E = 3 and compute D = E inv (mod B)  = B -
	//[B/3]
	ie = 3;
	id = ib - ib/3;
	if (idb >= 1) {
	    cout << "a = " << ia;
	    cout << "b = " << ib;
	    cout << "d = " << id;
	    cout << "e = " << ie;
	}
	time(&tm1);

	//!   Step 4: Encrypt and decrypt a series of numbers chosen
	//at random.
	for (int j = 1; j <= n1; j++) {
	    if (idb >= 1) cout << "Inner trial " << j << '\n';
	    ix = ia*mp_rand();
	    expm(idb, ix, ie, ia, iy);
	    expm(idb, iy, id, ia, iz);
	    if (ix != iz) {
		cout << "*** Test failed. " << j << '\n';
		cout << "x = " << ix;
		cout << "y = " << iy;
		cout << "z = " << iz;
		exit(1);
	    }
	}
	time(&tm2);
	if (idb >= 1) cout << "CPU times: " << difftime(tm1, tm0) <<
	  difftime(tm2, tm1) << '\n';
    }

    time(&tt1);
    cout << "Test Passed.\n";
//    cout << Total CPU time = " << difftime(tt1, tt0) << '\n';

    return 0;
}

