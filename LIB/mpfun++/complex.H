#ifndef	_MP_COMPLEX_
#define	_MP_COMPLEX_

#include <math.h>
#include <iostream.h>
#include "arch.h"

class FComplex;
class DComplex;

class FComplex {
    float re, im;
    const FComplex rtop() const
      {
	  return FComplex(FHYPOT(re, im), ATAN2F(im, re));
      }
    const FComplex ptor() const
      {
	  return FComplex(re*COSF(im), re*SINF(im));
      }
  public:
    // Constructors.
    FComplex()
      {
	  re = im = 0.0;
      }
    FComplex(float r, float i=0.0)
      {
	  re = r; im = i;
      }
    FComplex(const FComplex& c)
      {
	  re = c.re; im = c.im;
      }

    // Assignment.
    FComplex& operator=(const FComplex& c)
      {
	  re = c.re; im = c.im; return *this;
      }
    FComplex& operator=(float f)
      {
	  re = f; im = 0.0; return *this;
      }

    // Field extractors.
    friend const float real(const FComplex& c)
      {
	  return c.re;
      }
    friend const float imag(const FComplex& c)
      {
	  return c.im;
      }

    // Miscellaneous mathematical functions.
    friend const float abs(const FComplex& c)
      {
	  return FHYPOT(c.re, c.im);
      }
    friend const float arg(const FComplex& c)
      {
	  return ATAN2F(c.im, c.re);
      }
    friend const FComplex conj(const FComplex& c)
      {
	  return FComplex(c.re, -c.im);
      }
    friend const FComplex cos(const FComplex& c)
      {
	  // cos(x+iy) = cos(x)cosh(y) - i sin(x)sin(y) ??
	  cerr << "cos(const FComplex&): Not implemented\n";
	  return c;
      }
    friend const FComplex exp(const FComplex& c)
      {
	  float r = FEXP(c.re);
	  return FComplex(r*COSF(c.im), r*SINF(c.im));
      }
    friend const FComplex log(const FComplex& c)
      {
	  cerr << "log(const FComplex&): Not implemented\n";
	  return c;
      }
    friend const FComplex power(const float f, const FComplex& c)
      {
	  float r = POWF(f, c.re);
	  return FComplex(r*COSF(c.im), r*SINF(c.im));
      }
    friend const FComplex power(const FComplex& c, const int i)
      {
	  FComplex polar = c.rtop();
	  polar.re = POWF(polar.re, float(i));
	  polar.im *= i;
	  return polar.ptor();
      }
    friend const FComplex power(const FComplex& c, const float f)
      {
	  FComplex polar = c.rtop();
	  polar.re = POWF(polar.re, f);
	  polar.im *= f;
	  return polar.ptor();
      }
    friend const FComplex power(const FComplex& c1, const FComplex& c2)
      {
	  cerr << "power(const FComplex&, const FComplex&): Not implemented\n";
	  return c1;
      }
    friend const FComplex sin(const FComplex& c)
      {
	  // sin(x+iy) = sinx(x) cosh(y) + i cos(x) sinh(y)
	  cerr << "sin(const FComplex&): Not implemented\n";
	  return c;
      }
    friend const FComplex sqrt(const FComplex& c)
      {
	  FComplex polar = c.rtop();
	  polar.re = SQRTF(polar.re);
	  polar.im = 0.5*polar.im;
	  return polar.ptor();
      }

    // Overloaded arithmetic and logical operators.
    friend const FComplex operator+(const FComplex& c1, const FComplex& c2)
      {
	  return FComplex(c1.re+c2.re, c1.im+c2.im);
      }
    friend const FComplex operator-(const FComplex& c)
      {
	  return FComplex(-c.re, -c.im);
      }
    friend const FComplex operator-(const FComplex& c1, const FComplex& c2)
      {
	  return FComplex(c1.re-c2.re, c1.im-c2.im);
      }
    friend const FComplex operator*(const FComplex& c1, const FComplex& c2)
      {
	  return FComplex(c1.re*c2.re-c1.im*c2.im, c1.re*c2.im+c1.im*c2.re);
      }
    friend const FComplex operator/(const FComplex& c1,
				    const FComplex& c2)
      {
	  float denom = c2.re*c2.re+c2.im*c2.im;
	  return FComplex((c1.re*c2.re+c1.im*c2.im)/denom,
			  (c1.im*c2.re-c1.re*c2.im)/denom);
      }
    friend const int operator==(const FComplex& c1, const FComplex& c2)
      {
	  return (c1.re==c2.re && c1.im==c2.im) ? 1 : 0;
      }
    friend const int operator!=(const FComplex& c1, const FComplex& c2)
      {
	  return (c1.re!=c2.re || c1.im!=c2.im) ? 1 : 0;
      }

    // More overloaded opeators.
    void operator+=(const FComplex& c)
      {
	  re += c.re; im += c.im;
      }
    void operator-=(const FComplex& c)
      {
	  re -= c.re; im -= c.im;
      }
    void operator*=(const FComplex& c)
      {
	  re = re*c.re-im*c.im; im += im*c.re+re*c.im;
      }
    void operator/=(const FComplex& c)
      {
	  float denom = c.re*c.re+c.im*c.im;
	  re = (re*c.re+im*c.im)/denom;
	  im = (im*c.re-re*c.im)/denom;
      }

    // Stream I/O.
    friend istream& operator>>(istream& s, FComplex& c)
      {
	  s >> c.re >> c.im;
	  return s;
      }
    friend ostream& operator<<(ostream& s, const FComplex& c)
      {
	  s << '(' << c.re;
	  s << ',' << c.im ;
	  s << ")\n";
	  return s;
      }
};

class DComplex {
    double re, im;
    const DComplex rtop() const
      {
	  return DComplex(hypot(re, im), atan2(im, re));
      }
    const DComplex ptor() const
      {
	  return DComplex(re*cos(im), re*sin(im));
      }
  public:
    // Constructors.
    DComplex()
      {
	  re = im = 0.0;
      }
    DComplex(double r, double i=0.0)
      {
	  re = r; im = i;
      }
    DComplex(const DComplex& d)
      {
	  re = d.re; im = d.im;
      }
    DComplex(const FComplex& d)
      {
	  re = real(d); im = imag(d);
      }

    // Assignment.
    DComplex& operator=(const DComplex& d)
      {
	  re = d.re; im = d.im; return *this;
      }
    DComplex& operator=(double f)
      {
	  re = f; im = 0.0; return *this;
      }

    // Field extractors.
    friend const double real(const DComplex& d)
      {
	  return d.re;
      }
    friend const double imag(const DComplex& d)
      {
	  return d.im;
      }

    // Miscellaneous mathematical functions.
    friend const double abs(const DComplex& d)
      {
	  return hypot(d.re, d.im);
      }
    friend const double arg(const DComplex& d)
      {
	  return atan2(d.im, d.re);
      }
    friend const DComplex conj(const DComplex& d)
      {
	  return DComplex(d.re, -d.im);
      }
    friend const DComplex cos(const DComplex& d)
      {
	  cerr << "cos(const DComplex&): Not implemented\n";
	  return d;
      }
    friend const DComplex exp(const DComplex& d)
      {
	  float r = exp(d.re);
	  return DComplex(r*cos(d.im), r*sin(d.im));
      }
    friend const DComplex log(const DComplex& d)
      {
	  cerr << "log(const DComplex&): Not implemented\n";
	  return d;
      }
    friend const DComplex power(const double r, const DComplex& d)
      {
	  float t = pow(r, d.re);
	  return DComplex(t*cos(d.im), t*sin(d.im));
      }
    friend const DComplex power(const DComplex& d, const int i)
      {
	  DComplex polar = d.rtop();
	  polar.re = pow(polar.re, double(i));
	  polar.im *= i;
	  return polar.ptor();
      }
    friend const DComplex power(const DComplex& d, const double r)
      {
	  DComplex polar = d.rtop();
	  polar.re = pow(polar.re, r);
	  polar.im *= r;
	  return polar.ptor();
      }
    friend const DComplex power(const DComplex& d1, const DComplex& d2)
      {
	  cerr << "power(const DComplex&, const DComplex&): Not implemented\n";
	  return d1;
      }
    friend const DComplex sin(const DComplex& d)
      {
	  cerr << "sin(const DComplex&): Not implemented\n";
	  return d;
      }
    friend const DComplex sqrt(const DComplex& d)
      {
	  DComplex polar = d.rtop();
	  polar.re = sqrt(polar.re);
	  polar.im = 0.5*polar.im;
	  return polar.ptor();
      }

    // Overloaded arithmetic and logical operators.
    friend const DComplex operator+(const DComplex& d1, const DComplex& d2)
      {
	  return DComplex(d1.re+d2.re, d1.im+d2.im);
      }
    friend const DComplex operator-(const DComplex& d)
      {
	  return DComplex(-d.re, -d.im);
      }
    friend const DComplex operator-(const DComplex& d1, const DComplex& d2)
      {
	  return DComplex(d1.re-d2.re, d1.im-d2.im);
      }
    friend const DComplex operator*(const DComplex& d1, const DComplex& d2)
      {
	  return DComplex(d1.re*d2.re-d1.im*d2.im, d1.re*d2.im+d1.im*d2.re);
      }
    friend const DComplex operator/(const DComplex& d1,
				    const DComplex& d2)
      {
	  double denom = d2.re*d2.re+d2.im*d2.im;
	  return DComplex((d1.re*d2.re+d1.im*d2.im)/denom,
			  (d1.im*d2.re-d1.re*d2.im)/denom);
      }
    friend const int operator==(const DComplex& d1, const DComplex& d2)
      {
	  return (d1.re==d2.re && d1.im==d2.im) ? 1 : 0;
      }
    friend const int operator!=(const DComplex& d1, const DComplex& d2)
      {
	  return (d1.re!=d2.re || d1.im!=d2.im) ? 1 : 0;
      }

    // More overloaded operators.
    void operator+=(const DComplex& d)
      {
	  re += d.re; im += d.im;
      }
    void operator-=(const DComplex& d)
      {
	  re -= d.re; im -= d.im;
      }
    void operator*=(const DComplex& d)
      {
	  re = re*d.re-im*d.im; im += im*d.re+re*d.im;
      }
    void operator/=(const DComplex& d)
      {
	  double denom = d.re*d.re+d.im*d.im;
	  re = (re*d.re+im*d.im)/denom;
	  im = (im*d.re-re*d.im)/denom;
      }

    // Stream I/O.
    friend istream& operator>>(istream& s, DComplex& d)
      {
	  s >> d.re >> d.im;
	  return s;
      }
    friend ostream& operator<<(ostream& s, const DComplex& d)
      {
	  s << '(' << d.re;
	  s << ',' << d.im;
	  s << ")\n";
	  return s;
      }
};
#endif
