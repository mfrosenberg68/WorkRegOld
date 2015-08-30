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
// This package of C++ declarations provides a C++ interface to David H.
// Bailey's MPFUN multiprecision package written in standard Fortran-77.
// **************************************************************************

#include <iostream.h>
#include <fstream.h>
#include "complex.H"

const int mpipl = 100;
const int mpiou = 56;
const int mpiep = 10 - mpipl;

const int mpwds = ((mpipl * 100000) / 722472) + 1;
const int mp4 = mpwds + 4;
const int mp24 = 2 * mp4;
const int mp41 = mp4 + 1;

extern float	mpt0[mp24];
extern float	mpt1[mp24];
extern float	mpt2[mp24];
extern float	mpt3[mp24];
extern float	mpt4[mp24];

class mp {
  protected:
    void mpdexc(const char*, int, float*);
    void mpxzc(const DComplex&, float*);
    void mpmzc(const float*, float*);
  public:
    static char	az[mpipl+100];
    friend void mp_init();
};

class mp_integer;
class mp_real;
class mp_complex;

class mp_integer : public mp {
  public:
    float	mpi[mp4];

    mp_integer() {}
// c    mp_integer(mp_real&);
// c    mp_integer(mp_complex&);
    mp_integer(int);
    mp_integer(float);
//    mp_integer(complex);
    mp_integer(double);
//    mp_integer(dcomplex);
    mp_integer(char *);

    mp_integer& operator=(const int&);
    mp_integer& operator=(const float&);
    mp_integer& operator=(const double&);
    mp_integer& operator=(const char *);
    mp_integer& operator=(const mp_integer&);
    mp_integer& operator=(const mp_real&);
    mp_integer& operator=(const mp_complex&);

    mp_integer& operator+=(const mp_integer&);
    mp_integer& operator-=(const mp_integer&);
    mp_integer& operator*=(const mp_integer&);
    mp_integer& operator/=(const mp_integer&);

    friend const mp_integer operator+(const mp_integer&, const mp_integer&);
    friend const mp_integer operator-(const mp_integer&, const mp_integer&);
    friend const mp_integer operator-(const mp_integer&);
    friend const mp_integer operator*(const mp_integer&, const mp_integer&);
    friend const mp_integer operator/(const mp_integer&, const mp_integer&);
    friend const mp_integer power(const mp_integer&, const mp_integer&);
    friend const mp_real power(const mp_integer&, const mp_real&);
    friend const mp_integer power(const int, const mp_integer&);
    friend const mp_integer power(const mp_integer&, const int);
    friend const mp_real power(const float, const mp_integer&);
    friend const mp_real power(const mp_integer&, const float);
    friend const mp_real power(const double, const mp_integer&);
    friend const mp_real power(const mp_integer&, const double);
    friend int operator==(const mp_integer&, const mp_integer&);
    friend int operator!=(const mp_integer&, const mp_integer&);
    friend int operator<=(const mp_integer&, const mp_integer&);
    friend int operator>=(const mp_integer&, const mp_integer&);
    friend int operator<(const mp_integer&, const mp_integer&);
    friend int operator>(const mp_integer&, const mp_integer&);
    friend istream& operator>>(istream&, mp_integer&);
    friend ostream& operator<<(ostream&, const mp_integer&);
    friend ifstream& operator>>(ifstream&, mp_integer&);
    friend ofstream& operator<<(ofstream&, const mp_integer&);

    friend const mp_integer abs(const mp_integer&);
//    friend complex cmplx(const mp_integer&, const mp_integer&);
    friend double dble(const mp_integer&);
//    friend dcomplex dcmplx(const mp_integer&, const mp_integer&);
    friend int INT(const mp_integer&);
    friend const mp_integer max(const mp_integer&, const mp_integer&);
    friend const mp_integer min(const mp_integer&, const mp_integer&);
    friend const mp_integer mod(const mp_integer&, const mp_integer&);
    friend float real(const mp_integer&);
    friend const mp_integer sign(const mp_integer&, const mp_integer&);
};

class mp_real : public mp {
  public:
    float	mpr[mp4];

    mp_real() {}
    mp_real(const mp_integer&);
// c    mp_real(mp_complex&);
    mp_real(int);
    mp_real(float);
//    mp_real(complex);
    mp_real(double);
//    mp_real(dcomplex);
    mp_real(char *);

    mp_real& operator=(const int&);
    mp_real& operator=(const float&);
    mp_real& operator=(const double&);
    mp_real& operator=(const char *);
    mp_real& operator=(const mp_integer&);
    mp_real& operator=(const mp_real&);
    mp_real& operator=(const mp_complex&);

    mp_real& operator+=(const mp_real&);
    mp_real& operator-=(const mp_real&);
    mp_real& operator*=(const mp_real&);
    mp_real& operator/=(const mp_real&);

    friend const mp_real operator+(const mp_real&, const mp_real&);
    friend const mp_real operator-(const mp_real&, const mp_real&);
    friend const mp_real operator-(const mp_real&);
    friend const mp_real operator*(const mp_real&, const mp_real&);
    friend const mp_real operator/(const mp_real&, const mp_real&);
    friend const mp_real power(const mp_real&, const mp_integer&);
    friend const mp_real power(const mp_real&, const mp_real&);
    friend const mp_real power(const int, const mp_real&);
    friend const mp_real power(const mp_real&, const int);
    friend const mp_real power(const float, const mp_real&);
    friend const mp_real power(const mp_real&, const float);
    friend const mp_real power(const double, const mp_real&);
    friend const mp_real power(const mp_real&, const double);
    friend int operator==(const mp_real&, const mp_real&);
    friend int operator!=(const mp_real&, const mp_real&);
    friend int operator<=(const mp_real&, const mp_real&);
    friend int operator>=(const mp_real&, const mp_real&);
    friend int operator<(const mp_real&, const mp_real&);
    friend int operator>(const mp_real&, const mp_real&);
    friend istream& operator>>(istream&, mp_real&);
    friend ostream& operator<<(ostream&, const mp_real&);
    friend ifstream& operator>>(ifstream&, mp_real&);
    friend ofstream& operator<<(ofstream&, const mp_real&);

    friend const mp_real abs(const mp_real&);
    friend const mp_real acos(const mp_real&);
    friend const mp_real aint(const mp_real&);
    friend const mp_real anint(const mp_real&);
    friend const mp_real asin(const mp_real&);
    friend const mp_real atan(const mp_real&);
    friend const mp_real atan2(const mp_real&, const mp_real&);
//    friend complex cmplx(const mp_real&, const mp_real&);
    friend const mp_real cos(const mp_real&);
    friend const mp_real cosh(const mp_real&);
    friend double dble(const mp_real&);
//    friend dcomplex dcmplx(const mp_real&, const mp_real&);
    friend const mp_real exp(const mp_real&);
    friend int INT(const mp_real&);
    friend const mp_real log(const mp_real&);
    friend const mp_real log10(const mp_real&);
    friend const mp_real max(const mp_real&, const mp_real&);
    friend const mp_real min(const mp_real&, const mp_real&);
    friend const mp_real mod(const mp_real&, const mp_real&);
    friend void mpcsshf(const mp_real&, mp_real&, mp_real&);
    friend void mpcssnf(const mp_real&, mp_real&, mp_real&);
    friend const mp_real mpnrtf(const mp_real&, int);
    friend const mp_real mp_rand();
    friend const mp_integer nint(const mp_real&);
    friend float REAL(const mp_real&);
    friend const mp_real sign(const mp_real&, const mp_real&);
    friend const mp_real sin(const mp_real&);
    friend const mp_real sinh(const mp_real&);
    friend const mp_real sqrt(const mp_real&);
    friend const mp_real tan(const mp_real&);
    friend const mp_real tanh(const mp_real&);
};

class mp_complex : public mp {
  public:
    float	mpc[mp24];

    mp_complex() {}
    mp_complex(const mp_integer&);
    mp_complex(const mp_real&);
    mp_complex(int);
    mp_complex(float);
//    mp_complex(complex);
    mp_complex(double);
//    mp_complex(dcomplex);
    mp_complex(char *);
    mp_complex(const mp_integer&, const mp_integer&);
    mp_complex(const mp_real&, const mp_real&);
    mp_complex(int, int);
    mp_complex(float, float);
    mp_complex(double, double);
    mp_complex(char *, char *);

    mp_complex& operator=(const int&);
    mp_complex& operator=(const float&);
    mp_complex& operator=(const double&);
    mp_complex& operator=(const mp_integer&);
    mp_complex& operator=(const mp_real&);
    mp_complex& operator=(const mp_complex&);

    mp_complex& operator+=(const mp_complex&);
    mp_complex& operator-=(const mp_complex&);
    mp_complex& operator*=(const mp_complex&);
    mp_complex& operator/=(const mp_complex&);
    mp_complex& operator%=(const mp_complex&);

    friend const mp_complex operator+(const mp_complex&, const mp_complex&);
    friend const mp_complex operator-(const mp_complex&, const mp_complex&);
    friend const mp_complex operator-(const mp_complex&);
    friend const mp_complex operator*(const mp_complex&, const mp_complex&);
    friend const mp_complex operator/(const mp_complex&, const mp_complex&);
    friend const mp_complex power(const mp_complex&, const int);
    friend int operator==(const mp_complex&, const mp_complex&);
    friend int operator!=(const mp_complex&, const mp_complex&);
    friend istream& operator>>(istream&, mp_complex&);
    friend ostream& operator<<(ostream&, const mp_complex&);
    friend ifstream& operator>>(ifstream&, mp_complex&);
    friend ofstream& operator<<(ofstream&, const mp_complex&);

    friend const mp_real abs(const mp_complex&);
    friend const mp_real aimag(const mp_complex&);
//    friend complex cmplx(const mp_complex&);
    friend const mp_complex conjg(const mp_complex&);
    friend double dble(const mp_complex&);
//    friend dcomplex dcmplx(const mp_complex&);
    friend int INT(const mp_complex&);
    friend float REAL(const mp_complex&);
    friend mp_real MP_REAL(const mp_complex&);
    friend const mp_complex sqrt(const mp_complex&);
};

extern mp_real 	mpl02, mpl10, mppic, mpeps;
extern int	mpoud;

