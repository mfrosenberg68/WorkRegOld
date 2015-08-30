#include "mpmod.H"

mp_real dot(int n, mp_real* a, mp_real *b)

{
    mp_real	s = 0.0;

    for (int i = 1-1; i < n; i++)
	s = s + a[i]*b[i];
    return s;
}

int main()

/*
!  This is the test progrm for the F-90 based MP translation modules.
!
!  David H. Bailey   June 2, 1994
!
!  Translated into C++.
!
!  Siddhartha Chatterjee (RIACS) 
!  15 July 1994
*/

{
    mp_integer	ia, ib, ic;
    mp_complex	c, d, e;
    const int	n = 25;
    mp_real	a[n], b[n];

    mp_init();

    // Character-to-MP assignment, generic MPREAL function, pi and e.
    mp_real x = "1.234567890 1234567890 1234567890 D-100";
    mp_real ee = exp(mp_real(1.0));
    cout << x << mppic << ee;
    
    mp_real s = 0.0;
    // Loop with subscripted MP variables.
    for (int i = 1; i <= n; i++) {
	a[i-1] = 2*i+1;
	b[i-1] = 2.0*a[i-1]*(a[i-1]+1.0);
	s = s + power(b[i-1], 2);
    }
    cout << s;

    // Expression with mixed MPI and MPR entities.
    ia = s;
    ib = 262144;
    s = (s+327.25) * mod(ia, 4*ib);
    cout << s;

    // A complex square root reference.
    e = sqrt(mp_complex(2.0*s, s));
    cout << e;

    // External and intrinsic MP function references in expressions.
    s = dot(n, a, b);
    mp_real t = 2.0 * power(sqrt(s), 2);
    cout << s << t;

    s = s/1048576.0;
    t = s + 2.0*log(s);
    x = 3 + nint(t) * 5;
    cout << s << t << x;

    // A deeply nested expression with function references.
    x = (s + (2 * (s - 5) + 3 * (t - 5))) * exp(cos(log(s)));
    cout << x;

    // A special MP subroutine call (computes both cos and sin of S).
    mp_real	y;
    mpcssnf(s, x, y);
    t = 1.0 - (power(x, 2) + power(y, 2));

    // IF-THEN-ELSE construct involving MP variables.
    if (abs(t) < mpeps)
	cout << t;
    else
	cout << mpeps;
    return 0;
}
