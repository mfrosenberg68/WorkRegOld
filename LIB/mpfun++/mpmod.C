// **************************************************************************
// MPFUN++: A Multiple Precision Floating Point Computation Package
//
// C++ version
// Version Date: 22 February 1995
//
// Author:
//
//    Siddhartha Chatterjee		Telephone: +1 919 962 1766
//    Department of Computer Science	Facsimile: +1 919 962 1799
//    CB#3175, Sitterson Hall		Internet: sc@cs.unc.edu
//    The University of North Carolina
//    Chapel Hill, NC 27599-3175
//    USA
//
// Restrictions:
//
// Description:
//
// This package of C++ routines implements a C++ interface to David H.
// Bailey's MPFUN multiprecision package written in standarad Fortran-77.
// **************************************************************************

#include "mpmod.H"
#include "arch.h"

extern "C" {
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
}

//extern LANG float& SIGN(const float& a, const float& b);
extern LANG void MPADD(const float*, const float*, float*);
extern LANG void MPANG(const float*, const float*, const float*, float*);
extern LANG void MPCADD(const int&, const float*, const float*, float*);
extern LANG void MPCDIV(const int&, const float*, const float*, float*);
extern LANG void MPCEQ(const int&, const float*, float*);
extern LANG void MPCMUL(const int&, const float*, const float*, float*);
extern LANG void MPCPR(const float*, const float*, int&);
extern LANG void MPCPWR(const int&, const float*, const int&, float*);
extern LANG void MPCSQR(const int&, const float*, float*);
extern LANG void MPCSSH(const float*, const float*, float*, float*);
extern LANG void MPCSSN(const float*, const float*, float*, float*);
extern LANG void MPCSUB(const int&, const float*, const float*, float*);
extern LANG void MPDIV(const float*, const float*, float*);
extern LANG void MPDMC(const double&, const int&, float*);
extern LANG void MPEQ(const float*, float*);
extern LANG void MPEXP(const float*, const float*, float*);
extern LANG void MPINFR(const float*, float*, float*);
extern LANG void MPINPC(const char*, const int&, float*);
extern LANG void MPLOG(const float*, const float*, float*);
extern LANG void MPMDC(const float*, double&, int&);
extern LANG void MPMMPC(const float*, const float*, const int&, float*);
extern LANG void MPMUL(const float*, const float*, float*);
extern LANG void MPMULD(const float*, const double&, const int&, float*);
extern LANG void MPNINT(const float*, float*);
extern LANG void MPNPWR(const float*, const int&, float*);
extern LANG void MPNRT(const float*, const int&, float*);
extern LANG void MPOUT(const int&, const float*, const int&, char*);
extern LANG void MPOUTS(const float*, const int&, char*, int&);
extern LANG void MPPI(float*);
extern LANG void MPRAND(float*);
extern LANG void MPSETP(const char*, const int&, const int&);
extern LANG void MPSQRT(const float*, float*);
extern LANG void MPSUB(const float*, const float*, float*);

float	mpt0[mp24];
float	mpt1[mp24];
float	mpt2[mp24];
float	mpt3[mp24];
float	mpt4[mp24];

mp_real 	mpl02, mpl10, mppic, mpeps;
int		mpoud;
char 		mp::az[mpipl+100];

const float SIGN(const float& a, const float& b)

{
    return b < 0 ? -a : a;
}

void mp_init()

{
    MPSETP("nw", mpwds+1, 2);
    MPPI(mpt1);
    MPDMC(double(2), 0, mpt0);
    MPLOG(mpt0, mpt2, mpt2);
    MPDMC(double(10), 0, mpt0);
    MPLOG(mpt0, mpt2, mpt3);
    MPNPWR(mpt0, mpiep, mpt4);
    MPSETP("nw", mpwds, 2);
    MPEQ(mpt1, mppic.mpr);
    MPEQ(mpt2, mpl02.mpr);
    MPEQ(mpt3, mpl10.mpr);
    MPEQ(mpt4, mpeps.mpr);
    mpoud = mpiou;
}

void mp::mpdexc(const char *a, int l, float *b)

{
    for (int i = 0; i < l; i++) {
        if (a[i] == 'D' || a[i] == 'E' || a[i] == 'd' || a[i] == 'e')
          goto lab;
    }
    MPINPC(a, l, b);
    return;
  lab:
    int i1 = i;
    int l1 = i - 1;
    int l2 = l - i;
    char c[mpipl+100];
    c[0] = '1';
    c[1] = '0';
    c[2] = '^';
    char *cc = c+3;

    for (i = 1; i < l2; i++) {
	*cc = a[i+i1];
	cc++;
    }

    *cc = 'x'; cc++;

    for (i = 0; i <= l1; i++) {
	*cc = a[i];
	cc++;
    }

    *cc = '\0';

    int l12 = l1 + l2 + 4;
    MPINPC(c, l12, b);
}

void mp::mpxzc(const DComplex& a, float *b)

{
    double da = real(a);
    MPDMC(da, 0, b);
    da = imag(a);
    MPDMC(da, 0, b+mp41-1);
}

void mp::mpmzc(const float *a, float *b)

{
    MPEQ(a, b);
    b[mp41-1] = 0.0;
    b[mp4+2-1] = 1.0;
}

/***************** MP_INTEGER member functions **********************/

/*mp_integer::mp_integer(mp_real& qa)

{
    MPEQ((float*) qa.mpr, mpt1);
    MPINFR(mpt1, mpi, mpt2);
}
*/

/*mp_integer::mp_integer(mp_complex& za)

{
    MPEQ(za.mpc, mpt1);
    MPINFR(mpt1, mpi, mpt2);
}*/

mp_integer::mp_integer(int ia)

{
    MPDMC(double(ia), 0, mpi);
}

mp_integer::mp_integer(float ra)

{
    MPDMC(double(ra), 0, mpt1);
    MPINFR(mpt1, mpi, mpt2);
}

mp_integer::mp_integer(double da)

{
    MPDMC(da, 0, mpt1);
    MPINFR(mpt1, mpi, mpt2);
}

mp_integer& mp_integer::operator=(const int& ib)

{
    MPDMC(double(ib), 0, mpi);
    return *this;
}

mp_integer& mp_integer::operator=(const float& rb)

{
    MPDMC(double(rb), 0, mpt1);
    MPINFR(mpt1, mpi, mpt2);
    return *this;
}

mp_integer& mp_integer::operator=(const double& db)

{
    MPDMC(db, 0, mpt1);
    MPINFR(mpt1, mpi, mpt2);
    return *this;
}

mp_integer& mp_integer::operator=(const mp_integer& jb)

{
    MPEQ(jb.mpi, mpi);
    return *this;
}

mp_integer& mp_integer::operator=(const mp_real& qb)

{
    MPEQ(qb.mpr, mpt1);
    MPINFR(mpt1, mpi, mpt2);
    return *this;
}

mp_integer& mp_integer::operator=(const mp_complex& zb)

{
    MPEQ(zb.mpc, mpt1);
    MPINFR(mpt1, mpi, mpt2);
    return *this;
}

mp_integer& mp_integer::operator=(const char * ab)

{
    mpdexc(ab, strlen(ab), mpt1);
    MPINFR(mpt1, mpi, mpt2);
    return *this;
}

mp_integer& mp_integer::operator+=(const mp_integer& ja)

{
    *this = (*this)+ja;
    return *this;
}

mp_integer& mp_integer::operator-=(const mp_integer& ja)

{
    *this = (*this)-ja;
    return *this;
}

mp_integer& mp_integer::operator*=(const mp_integer& ja)

{
    *this = (*this)*ja;
    return *this;
}

mp_integer& mp_integer::operator/=(const mp_integer& ja)

{
    *this = (*this)/ja;
    return *this;
}

const mp_integer operator+(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer res;
    
    MPADD(ja.mpi, jb.mpi, mpt1);
    MPINFR(mpt1, res.mpi, mpt2);
    return res;
}

const mp_integer operator-(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer res;

    MPSUB(ja.mpi, jb.mpi, mpt1);
    MPINFR(mpt1, res.mpi, mpt2);
    return res;
}

const mp_integer operator-(const mp_integer& ja)

{
    mp_integer res;

    MPEQ(ja.mpi, res.mpi);
    res.mpi[1-1] = - ja.mpi[1-1];
    return res;
}

const mp_integer operator*(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer res;
    
    MPMUL(ja.mpi, jb.mpi, mpt1);
    MPINFR(mpt1, res.mpi, mpt2);
    return res;
}

const mp_integer operator/(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer res;
    
    MPDIV(ja.mpi, jb.mpi, mpt1);
    MPINFR(mpt1, res.mpi, mpt2);
    return res;
}

const mp_integer power(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer res;

    MPLOG(ja.mpi, mpl02.mpr, mpt1);
    MPMUL(mpt1, jb.mpi, mpt2);
    MPEXP(mpt2, mpl02.mpr, mpt1);
    MPNINT(mpt1, res.mpi);
    return res;
}

const mp_real power(const mp_integer& ja, const mp_real& qb)

{
    mp_real res;
    
    MPLOG(ja.mpi, mpl02.mpr, mpt1);
    MPMUL(mpt1, qb.mpr, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

const mp_integer power(const int ia, const mp_integer& jb)

{
    mp_integer res;

    MPDMC(double(ia), 0, mpt1);
    MPLOG(mpt1, mpl02.mpr, mpt2);
    MPMUL(mpt2, jb.mpi, mpt3);
    MPEXP(mpt3, mpl02.mpr, mpt1);
    MPNINT(mpt1, res.mpi);
    return res;
}

const mp_integer power(const mp_integer& ja, const int ib)

{
    mp_integer res;
    
    MPNPWR(ja.mpi, ib, mpt1);
    MPNINT(mpt1, res.mpi);
    return res;
}

const mp_real power(const float ra, const mp_integer& jb)

{
    mp_real res;

    MPDMC(double(ra), 0, mpt1);
    MPLOG(mpt1, mpl02.mpr, mpt2);
    MPMUL(mpt2, jb.mpi, mpt3);
    MPEXP(mpt3, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_integer& ja, const float rb)

{
    mp_real res;

    MPLOG(ja.mpi, mpl02.mpr, mpt1);
    MPMULD(mpt1, double(rb), 0, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const double da, const mp_integer& jb)

{
    mp_real res;
    
    MPDMC(da, 0, mpt1);
    MPLOG(mpt1, mpl02.mpr, mpt2);
    MPMUL(mpt2, jb.mpi, mpt3);
    MPEXP(mpt3, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_integer& ja, const double db)

{
    mp_real res;

    MPLOG(ja.mpi, mpl02.mpr, mpt1);
    MPMULD(mpt1, db, 0, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

int operator==(const mp_integer& ja, const mp_integer& jb)

{
    int ic;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    return ic == 0 ? 1 : 0;
}

int operator!=(const mp_integer& ja, const mp_integer& jb)

{
    int ic;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    return ic != 0 ? 1 : 0;
}

int operator<=(const mp_integer& ja, const mp_integer& jb)

{
    int ic;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    return ic <= 0 ? 1 : 0;
}

int operator>=(const mp_integer& ja, const mp_integer& jb)

{
    int ic;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    return ic >= 0 ? 1 : 0;
}

int operator<(const mp_integer& ja, const mp_integer& jb)

{
    int ic;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    return ic < 0 ? 1 : 0;
}

int operator>(const mp_integer& ja, const mp_integer& jb)

{
    int ic;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    return ic > 0 ? 1 : 0;
}

istream& operator>>(istream& s, mp_integer& ja)

{
    int i = 0;
    char c = ' ';
    while (1) {
	if (! s >> c) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    ja.mpdexc(mp::az, i, ja.mpi);
    return s;
}

ostream& operator<<(ostream& s, const mp_integer& ja)

{
    int	nd = 0;
    for (int i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(ja.mpi, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s << mp::az[i];
    s << '\n';
    s.flush();
    return s;
}

ifstream& operator>>(ifstream& s, mp_integer& ja)

{
    int i = 0;
    char c;
    while (1) {
	if (! s.get(c)) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    ja.mpdexc(mp::az, i, ja.mpi);
    return s;
}

ofstream& operator<<(ofstream& s, mp_integer& ja)

{
    int	nd = 0;
    for (int i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(ja.mpi, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s.put(mp::az[i]);
    s.put('\n');
    s.flush();
    return s;
}

const mp_integer abs(const mp_integer& ja)

{
    mp_integer res;
    
    MPEQ(ja.mpi, res.mpi);
    res.mpi[1-1] = fabs(double(ja.mpi[1-1]));
    return res;
}

//complex cmplx(const mp_integer&, const mp_integer&);

double dble(const mp_integer& ja) 

{
    double 	da;
    int		ia;
    
    MPMDC(ja.mpi, da, ia);
    return da * pow(double(2), ia);
}

//dcomplex mp_integer::dcmplx(const mp_integer&, const mp_integer&);

int INT(const mp_integer& ja)

{
    double 	da;
    int		ia;
    
    MPMDC(ja.mpi, da, ia);
    return int(da * pow(double(2), ia));
}

const mp_integer max(const mp_integer& ja, const mp_integer& jb)

{
    int		ic;
    mp_integer res;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    if (ic >= 0)
        MPEQ(ja.mpi, res.mpi);
    else
        MPEQ(jb.mpi, res.mpi);
    return res;
}

const mp_integer min(const mp_integer& ja, const mp_integer& jb)

{
    int		ic;
    mp_integer 	res;
    
    MPCPR(ja.mpi, jb.mpi, ic);
    if (ic < 0)
        MPEQ(ja.mpi, res.mpi);
    else
        MPEQ(jb.mpi, res.mpi);
    return res;
}

const mp_integer mod(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer	res;

    MPDIV(ja.mpi, jb.mpi, mpt1);
    MPINFR(mpt1, mpt2, mpt3);
    MPMUL(jb.mpi, mpt2, mpt1);
    MPSUB(ja.mpi, mpt1, res.mpi);
    return res;
}

float real(const mp_integer& ja)

{
    double	da;
    int		ia;
    
    MPMDC(ja.mpi, da, ia);
    return da * pow(double(2), ia);
}

const mp_integer sign(const mp_integer& ja, const mp_integer& jb)

{
    mp_integer res;
    
    MPEQ(ja.mpi, res.mpi);
    res.mpi[1-1] = SIGN(res.mpi[1-1], jb.mpi[1-1]);
    return res;
}

/***************** MP_REAL member functions **********************/

mp_real::mp_real(const mp_integer& ja)

{
    MPEQ(ja.mpi, mpr);
}

/*mp_real::mp_real(mp_complex& za)

{
    MPEQ(za.mpc, mpr);
}*/

mp_real::mp_real(int ia)

{
    MPDMC(double(ia), 0, mpr);
}

mp_real::mp_real(float ra)

{
    MPDMC(double(ra), 0, mpr);
}

//mp_real(complex)

mp_real::mp_real(double da)

{
    MPDMC(da, 0, mpr);
}

//mp_real(dcomplex)

mp_real::mp_real(char * aa)

{
    int l = strlen(aa);
    
    for (int i = 0; i < l; i++) az[i] = aa[i];
    mpdexc(az, l, mpr);
}

mp_real& mp_real::operator=(const int& ib)

{
    MPDMC(double(ib), 0, mpr);
    return *this;
}

mp_real& mp_real::operator=(const float& rb)

{
    MPDMC(double(rb), 0, mpr);
    return *this;
}

mp_real& mp_real::operator=(const double& db)

{
    MPDMC(db, 0, mpr);
    return *this;
}

mp_real& mp_real::operator=(const mp_integer& jb)

{
    MPEQ(jb.mpi, mpr);
    return *this;
}

mp_real& mp_real::operator=(const mp_real& qb)

{
    MPEQ(qb.mpr, mpr);
    return *this;
}

mp_real& mp_real::operator=(const mp_complex& zb)

{
    MPEQ(zb.mpc, mpr);
    return *this;
}

mp_real& mp_real::operator=(const char * ab)

{
    mpdexc(ab, strlen(ab), mpr);
    return *this;
}

mp_real& mp_real::operator+=(const mp_real& qa)

{
    *this = (*this)+qa;
    return *this;
}

mp_real& mp_real::operator-=(const mp_real& qa)

{
    *this = (*this)-qa;
    return *this;
}

mp_real& mp_real::operator*=(const mp_real& qa)

{
    *this = (*this)*qa;
    return *this;
}

mp_real& mp_real::operator/=(const mp_real& qa)

{
    *this = (*this)/qa;
    return *this;
}

const mp_real operator+(const mp_real& qa, const mp_real& qb)

{
    mp_real res;
    
    MPADD(qa.mpr, qb.mpr, res.mpr);
    return res;
}

const mp_real operator-(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;

    MPSUB(qa.mpr, qb.mpr, res.mpr);
    return res;
}

const mp_real operator-(const mp_real& qa)

{
    mp_real	res;

    MPEQ(qa.mpr, res.mpr);
    res.mpr[1-1] = - qa.mpr[1-1];
    return res;
}

const mp_real operator*(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;

    MPMUL(qa.mpr, qb.mpr, res.mpr);
    return res;
}

const mp_real operator/(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;

    MPDIV(qa.mpr, qb.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_real& qa, const mp_integer& jb)

{
    mp_real res;

    MPLOG(qa.mpr, mpl02.mpr, mpt1);
    MPMUL(mpt1, jb.mpi, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;
    
    MPLOG(qa.mpr, mpl02.mpr, mpt1);
    MPMUL(mpt1, qb.mpr, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const int ia, const mp_real& qb)

{
    mp_real res;
    
    MPDMC(double(ia), 0, mpt1);
    MPLOG(mpt1, mpl02.mpr, mpt2);
    MPMUL(mpt2, qb.mpr, mpt3);
    MPEXP(mpt3, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_real& qa, const int ib)

{
    mp_real res;

    MPNPWR(qa.mpr, ib, res.mpr);
    return res;
}

const mp_real power(const float ra, const mp_real& qb)

{
    mp_real res;

    MPDMC(double(ra), 0, mpt1);
    MPLOG(mpt1, mpl02.mpr, mpt2);
    MPMUL(mpt2, qb.mpr, mpt3);
    MPEXP(mpt3, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_real& qa, const float rb)

{
    mp_real res;
    
    MPLOG(qa.mpr, mpl02.mpr, mpt1);
    MPMULD(mpt1, double(rb), 0, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const double da, const mp_real& qb)

{
    mp_real res;
    
    MPDMC(da, 0, mpt1);
    MPLOG(mpt1, mpl02.mpr, mpt2);
    MPMUL(mpt2, qb.mpr, mpt3);
    MPEXP(mpt3, mpl02.mpr, res.mpr);
    return res;
}

const mp_real power(const mp_real& qa, const double db)

{
    mp_real res;

    MPLOG(qa.mpr, mpl02.mpr, mpt1);
    MPMULD(mpt1, db, 0, mpt2);
    MPEXP(mpt2, mpl02.mpr, res.mpr);
    return res;
}

int operator==(const mp_real& qa, const mp_real& qb)

{
    int ic;
    
    MPCPR(qa.mpr, qb.mpr, ic);
    return ic == 0 ? 1 : 0;
}

int operator!=(const mp_real& qa, const mp_real& qb)

{
    int ic;
    
    MPCPR(qa.mpr, qb.mpr, ic);
    return ic != 0 ? 1 : 0;
}

int operator<=(const mp_real& qa, const mp_real& qb)

{
    int ic;
    
    MPCPR(qa.mpr, qb.mpr, ic);
    return ic <= 0 ? 1 : 0;
}

int operator>=(const mp_real& qa, const mp_real& qb)

{
    int ic;
    
    MPCPR(qa.mpr, qb.mpr, ic);
    return ic >= 0 ? 1 : 0;
}

int operator<(const mp_real& qa, const mp_real& qb)

{
    int ic;
    
    MPCPR(qa.mpr, qb.mpr, ic);
    return ic < 0 ? 1 : 0;
}

int operator>(const mp_real& qa, const mp_real& qb)

{
    int ic;
    
    MPCPR(qa.mpr, qb.mpr, ic);
    return ic > 0 ? 1 : 0;
}

istream& operator>>(istream& s, mp_real& qa)

{
    int i = 0;
    char c = ' ';
    while (1) {
	if (! s >> c) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    qa.mpdexc(mp::az, i, qa.mpr);
    return s;
}

ostream& operator<<(ostream& s, const mp_real& qa)

{
    int	nd = 0;
    for (int i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(qa.mpr, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s << mp::az[i];
    s << '\n';
    s.flush();
    return s;
}

ifstream& operator>>(ifstream& s, mp_real& qa)

{
    int i = 0;
    char c;
    while (1) {
	if (! s.get(c)) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    qa.mpdexc(mp::az, i, qa.mpr);
    return s;
}

ofstream& operator<<(ofstream& s, const mp_real& qa)

{
    int	nd = 0;
    for (int i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(qa.mpr, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s.put(mp::az[i]);
    s.put('\n');
    s.flush();
    return s;
}

const mp_real abs(const mp_real& qa)

{
    mp_real	res;

    MPEQ(qa.mpr, res.mpr);
    res.mpr[1-1] = fabs(double(qa.mpr[1-1]));
    return res;
}

const mp_real acos(const mp_real& qa)

{
    mp_real	res;

    MPDMC(double(1), 0, mpt1);
    MPMUL(qa.mpr, qa.mpr, mpt2);
    MPSUB(mpt1, mpt2, mpt3);
    MPSQRT(mpt3, mpt1);
    MPANG(qa.mpr, mpt1, mppic.mpr, res.mpr);
    return res;
}

const mp_real aint(const mp_real& qa)

{
    mp_real	res;

    MPINFR(qa.mpr, res.mpr, mpt1);
    return res;
}

const mp_real anint(const mp_real& qa)

{
    mp_real	res;

    MPNINT(qa.mpr, res.mpr);
    return res;
}

const mp_real asin(const mp_real& qa)

{
    mp_real	res;
    
    MPDMC(double(1), 0, mpt1);
    MPMUL(qa.mpr, qa.mpr, mpt2);
    MPSUB(mpt1, mpt2, mpt3);
    MPSQRT(mpt3, mpt1);
    MPANG(mpt1, qa.mpr, mppic.mpr, res.mpr);
    return res;
}

const mp_real atan(const mp_real& qa)

{
    mp_real	res;

    MPDMC(double(1), 0, mpt1);
    MPANG(qa.mpr, mpt1, mppic.mpr, res.mpr);
    return res;
}

const mp_real atan2(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;

    MPANG(qb.mpr, qa.mpr, mppic.mpr, res.mpr);
    return res;
}

//complex cmplx(const mp_real& qa, const mp_real& qb)

const mp_real cos(const mp_real& qa)

{
    mp_real	res;

    MPCSSN(qa.mpr, mppic.mpr, res.mpr, mpt1);
    return res;
}

const mp_real cosh(const mp_real& qa)

{
    mp_real	res;

    MPCSSH(qa.mpr, mpl02.mpr, res.mpr, mpt1);
    return res;
}

double dble(const mp_real& qa)

{
    double	da;
    int		ia;
    
    MPMDC(qa.mpr, da, ia);
    return da * pow(double(2), ia);
}

//dcomplex dcmplx(const mp_real& qa, const mp_real& qb)

const mp_real exp(const mp_real& qa)

{
    mp_real	res;

    MPEXP(qa.mpr, mpl02.mpr, res.mpr);
    return res;
}

int INT(const mp_real& qa)

{
    double	da;
    int		ia;

    MPMDC(qa.mpr, da, ia);
    return int(da * pow(double(2), ia));
}

const mp_real log(const mp_real& qa)

{
    mp_real	res;

    MPLOG(qa.mpr, mpl02.mpr, res.mpr);
    return res;
}

const mp_real log10(const mp_real& qa)

{
    mp_real	res;

    MPLOG(qa.mpr, mpl02.mpr, res.mpr);
    MPDIV(mpt1, mpl10.mpr, res.mpr);
    return res;
}

const mp_real max(const mp_real& qa, const mp_real& qb)

{
    int		ic;
    mp_real	res;

    MPCPR(qa.mpr, qb.mpr, ic);
    if (ic >= 0)
        MPEQ(qa.mpr, res.mpr);
    else
        MPEQ(qb.mpr, res.mpr);
    return res;
}

const mp_real min(const mp_real& qa, const mp_real& qb)

{
    int		ic;
    mp_real	res;

    MPCPR(qa.mpr, qb.mpr, ic);
    if (ic < 0)
        MPEQ(qa.mpr, res.mpr);
    else
        MPEQ(qb.mpr, res.mpr);
    return res;
}

const mp_real mod(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;
    
    MPDIV(qa.mpr, qb.mpr, mpt1);
    MPINFR(mpt1, mpt2, mpt3);
    MPMUL(qb.mpr, mpt2, mpt1);
    MPSUB(qa.mpr, mpt1, res.mpr);
    return res;
}

void mpcsshf(const mp_real& qa, mp_real& qb, mp_real& qc)

{
    MPCSSH(qa.mpr, mpl02.mpr, qb.mpr, qc.mpr);
}

void mpcssnf(const mp_real& qa, mp_real& qb, mp_real& qc)

{
    MPCSSN(qa.mpr, mppic.mpr, qb.mpr, qc.mpr);
}

const mp_real mpnrtf(const mp_real& qa, int ib)

{
    mp_real	res;
    
    MPNRT(qa.mpr, ib, res.mpr);
    return res;
}

const mp_real mp_rand()

{
    mp_real	res;

    MPRAND(res.mpr);
    return res;
}

const mp_integer nint(const mp_real& qa)

{
    mp_integer	res;

    MPNINT(qa.mpr, res.mpi);
    return res;
}

float REAL(const mp_real& qa)

{
    double	da;
    int		ia;

    MPMDC(qa.mpr, da, ia);
    return da * pow(double(2), ia);
}

const mp_real sign(const mp_real& qa, const mp_real& qb)

{
    mp_real	res;

    MPEQ(qa.mpr, res.mpr);
    res.mpr[1-1] = SIGN(res.mpr[1-1], qb.mpr[1-1]);
    return res;
}

const mp_real sin(const mp_real& qa)

{
    mp_real	res;

    MPCSSN(qa.mpr, mppic.mpr, mpt1, res.mpr);
    return res;
}

const mp_real sinh(const mp_real& qa)

{
    mp_real	res;

    MPCSSH(qa.mpr, mpl02.mpr, mpt1, res.mpr);
    return res;
}

const mp_real sqrt(const mp_real& qa)

{
    mp_real	res;

    MPSQRT(qa.mpr, res.mpr);
    return res;
}

const mp_real tan(const mp_real& qa)

{
    mp_real	res;

    MPCSSN(qa.mpr, mppic.mpr, mpt1, mpt2);
    MPDIV(mpt1, mpt2, res.mpr);
    return res;
}

const mp_real tanh(const mp_real& qa)

{
    mp_real	res;

    MPCSSH(qa.mpr, mppic.mpr, mpt1, mpt2);
    MPDIV(mpt1, mpt2, res.mpr);
    return res;
}

/***************** MP_COMPLEX member functions **********************/

mp_complex::mp_complex(const mp_integer& ja)

{
    mpmzc(ja.mpi, mpc);
}

mp_complex::mp_complex(const mp_real& qa)

{
    mpmzc(qa.mpr, mpc);
}

mp_complex::mp_complex(int ia)

{
    mpxzc(DComplex(ia), mpc);
}

mp_complex::mp_complex(float ra)

{
    mpxzc(DComplex(ra), mpc);
}

//mp_complex::mp_complex(complex ca)

mp_complex::mp_complex(double da)

{
    mpxzc(DComplex(da), mpc);
}

//mp_complex::mp_complex(DComplex xa)

mp_complex::mp_complex(char* aa)

{
    int	l = strlen(aa);
    
    for (int i = 0; i < l; i++) {
        az[i] = aa[i];
    }
    MPINPC(az, l, mpt1);
    mpmzc(mpt1, mpc);
}

mp_complex::mp_complex(const mp_integer& ja, const mp_integer& jb)

{
    MPMMPC(ja.mpi, jb.mpi, mp4, mpc);
}

mp_complex::mp_complex(const mp_real& qa, const mp_real& qb)

{
    MPMMPC(qa.mpr, qb.mpr, mp4, mpc);
}

mp_complex::mp_complex(int ia, int ib)

{
    DComplex xa(ia, ib);
    mpxzc(xa, mpc);
}

mp_complex::mp_complex(float ra, float rb)

{
    DComplex xa(ra, rb);
    mpxzc(xa, mpc);
}

mp_complex::mp_complex(double da, double db)

{
    DComplex xa(da, db);
    mpxzc(xa, mpc);
}

mp_complex::mp_complex(char* aa, char* ab)

{
    int	l = strlen(aa);

    for (int i = 0; i < l; i++) {
        az[i] = aa[i];
    }
    MPINPC(az, l, mpc);
    l = strlen(ab);
    for (i = 0; i < l; i++) {
        az[i] = ab[i];
    }
    MPINPC(az, l, mpc+mp41-1);
}

mp_complex& mp_complex::operator=(const int& ib)

{
    DComplex xb = ib;
    mpxzc(xb, mpc);
    return *this;
}

mp_complex& mp_complex::operator=(const float& rb)

{
    DComplex xb = rb;
    mpxzc(xb, mpc);
    return *this;
}

mp_complex& mp_complex::operator=(const double& db)

{
    DComplex xb = db;
    mpxzc(xb, mpc);
    return *this;
}

mp_complex& mp_complex::operator=(const mp_integer& jb)

{
    mpmzc(jb.mpi, mpc);
    return *this;
}

mp_complex& mp_complex::operator=(const mp_real& qb)

{
    mpmzc(qb.mpr, mpc);
    return *this;
}

mp_complex& mp_complex::operator=(const mp_complex& zb)
{
    MPCEQ(mp4, zb.mpc, mpc);
    return *this;
}

mp_complex& mp_complex::operator+=(const mp_complex& za)

{
    *this = (*this)+za;
    return *this;
}

mp_complex& mp_complex::operator-=(const mp_complex& za)

{
    *this = (*this)-za;
    return *this;
}

mp_complex& mp_complex::operator*=(const mp_complex& za)

{
    *this = (*this)*za;
    return *this;
}

mp_complex& mp_complex::operator/=(const mp_complex& za)

{
    *this = (*this)/za;
    return *this;
}

const mp_complex operator+(const mp_complex& za, const mp_complex& zb)

{
    mp_complex	res;
    
    MPCADD(mp4, za.mpc, zb.mpc, res.mpc);
    return res;
}

const mp_complex operator-(const mp_complex& za, const mp_complex& zb)

{
    mp_complex	res;
    
    MPCSUB(mp4, za.mpc, zb.mpc, res.mpc);
    return res;
}

const mp_complex operator-(const mp_complex& za)

{
    mp_complex	res;
    
    MPCEQ(mp4, za.mpc, res.mpc);
    res.mpc[1-1] = - za.mpc[1-1];
    res.mpc[mp41-1] = - za.mpc[mp41-1];
    return res;
}

const mp_complex operator*(const mp_complex& za, const mp_complex& zb)

{
    mp_complex	res;
    
    MPCMUL(mp4, za.mpc, zb.mpc, res.mpc);
    return res;
}

const mp_complex operator/(const mp_complex& za, const mp_complex& zb)

{
    mp_complex	res;
    
    MPCDIV(mp4, za.mpc, zb.mpc, res.mpc);
    return res;
}

const mp_complex power(const mp_complex& za, const int ib)

{
    mp_complex	res;

    MPCPWR(mp4, za.mpc, ib, res.mpc);
    return res;
}

int operator==(const mp_complex& za, const mp_complex& zb)

{
    int	ic1, ic2;
    
    MPCPR(za.mpc, zb.mpc, ic1);
    MPCPR(za.mpc + mp41 - 1, zb.mpc + mp41 - 1, ic2);
    return (ic1 == 0 && ic2 == 0) ? 1 : 0;
}

int operator!=(const mp_complex& za, const mp_complex& zb)

{
    int	ic1, ic2;
    
    MPCPR(za.mpc, zb.mpc, ic1);
    MPCPR(za.mpc + mp41 - 1, zb.mpc + mp41 - 1, ic2);
    return (ic1 != 0 || ic2 != 0) ? 1 : 0;
}

istream& operator>>(istream& s, mp_complex& za)

{
    int i = 0;
    char c = ' ';
    while (1) {
	if (! s >> c) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    za.mpdexc(mp::az, i, za.mpc);
    while (1) {
	if (! s >> c) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    za.mpdexc(mp::az, i, za.mpc+mp41-1);
    return s;
}

ostream& operator<<(ostream& s, const mp_complex& za)

{
    int nd = 0;
    for (int i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(za.mpc, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s << mp::az[i];
    s << '\n';
    for (i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(za.mpc+mp41-1, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s << mp::az[i];
    s << '\n';
    s.flush();
    return s;
}

ifstream& operator>>(ifstream& s, mp_complex& za)

{
    int i = 0;
    char c;
    while (1) {
	if (! s.get(c)) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    za.mpdexc(mp::az, i, za.mpc);
    while (1) {
	if (! s.get(c)) return s;
	if (c == ',') break;
	if (c == '\n' || c == '\t' || c == ' ') continue;
	mp::az[i++] = c;
    }
    za.mpdexc(mp::az, i, za.mpc+mp41-1);
    return s;
}

ofstream& operator<<(ofstream& s, mp_complex& za)

{
    int nd = 0;
    for (int i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(za.mpc, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s.put(mp::az[i]);
    s.put('\n');
    for (i = 0; i < mpipl+100; i++) mp::az[i] = '$';
    MPOUTS(za.mpc+mp41-1, mpoud, mp::az, nd);
    assert(mp::az[nd-1] == ',');
    for (i = 0; i < nd; i++) s.put(mp::az[i]);
    s.put('\n');
    s.flush();
    return s;
}

const mp_real abs(const mp_complex& za)

{
    mp_real res;

    MPMUL(za.mpc, za.mpc, mpt1);
    MPMUL(za.mpc+mp41-1, za.mpc+mp41-1, mpt2);
    MPADD(mpt1, mpt2, mpt3);
    MPSQRT(mpt3, res.mpr);
    return res;
}

const mp_real aimag(const mp_complex& za)

{
    mp_real	res;
    
    MPEQ(za.mpc+mp41-1, res.mpr);
    return res;
}

//complex cmplx(const mp_complex& za)

const mp_complex conjg(const mp_complex& za)

{
    mp_complex	res;

    MPCEQ(mp4, za.mpc, res.mpc);
    res.mpc[mp41-1] = - za.mpc[mp41-1];
    return res;
}

double dble(const mp_complex& za)

{
    double	da;
    int		ia;

    MPMDC(za.mpc, da, ia);
    return da * pow(double(2), ia);
}

//DComplex dcmplx(const mp_complex& za)

int INT(const mp_complex& za)

{
    double	da;
    int		ia;

    MPMDC(za.mpc, da, ia);
    return int(da * pow(double(2), ia));
}

float REAL(const mp_complex& za)

{
    double	da;
    int		ia;

    MPMDC(za.mpc, da, ia);
    return da * pow(double(2), ia);
}

mp_real MP_REAL(const mp_complex& za)

{
    mp_real res;
    MPEQ(za.mpc, res.mpr);
    return res;
}

const mp_complex sqrt(const mp_complex& za)

{
    mp_complex	res;
    
    MPCSQR(mp4, za.mpc, res.mpc);
    return res;
}



