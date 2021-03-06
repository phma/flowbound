/******************************************************/
/*                                                    */
/* eisenstein.h - Eisenstein integers                 */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#ifndef EISENSTEIN_H
#define EISENSTEIN_H

#include <cmath>
#include <cstdlib>
#include <map>
#include <complex>
#include <gmpxx.h>
#define M_SQRT_3_4 0.86602540378443864676372317
// The continued fraction expansion is 0;1,6,2,6,2,6,2,...
#define M_SQRT_3 1.73205080756887729352744634
#define M_SQRT_1_3 0.5773502691896257645091487805

//PAGERAD should be 1 or 6 mod 8, which makes PAGESIZE 7 mod 8.
// The maximum is 147 because of the file format.
#define PAGERAD 6
#define PAGESIZE (PAGERAD*(PAGERAD+1)*3+1)
#define sqr(a) ((a)*(a))
extern const std::complex<double> omega,ZLETTERMOD;

class Eisenstein
{
private:
  mpz_class x,y; // x is the real part, y is at 120°
  static mpz_class numx,numy,denx,deny,quox,quoy,remx,remy;
  void divmod(Eisenstein b);
public:
  Eisenstein()
  {
    x=y=0;
  }
  Eisenstein(int xa)
  {
    x=xa;
    y=0;
  }
  Eisenstein(int xa,int ya)
  {
    x=xa;
    y=ya;
  }
  Eisenstein(mpz_class xa)
  {
    x=xa;
    y=0;
  }
  Eisenstein(mpz_class xa,mpz_class ya)
  {
    x=xa;
    y=ya;
  }
  Eisenstein(std::complex<double> z);
  Eisenstein operator+(Eisenstein b);
  Eisenstein operator-();
  Eisenstein operator-(Eisenstein b);
  Eisenstein operator*(Eisenstein b);
  Eisenstein& operator*=(Eisenstein b);
  Eisenstein& operator+=(Eisenstein b);
  Eisenstein operator/(Eisenstein b);
  Eisenstein operator%(Eisenstein b);
  bool operator==(Eisenstein b);
  bool operator!=(Eisenstein b);
  friend bool operator<(const Eisenstein a,const Eisenstein b); // only for the map
  mpz_class norm();
  int pageinx(int size,int nelts);
  int pageinx();
  int letterinx();
  mpz_class getx()
  {
    return x;
  }
  mpz_class gety()
  {
    return y;
  }
  mpz_class cartx()
  { // This returns twice the actual x coordinate.
    return 2*x-y;
  }
  mpz_class carty()
  { // This returns √(4/3) times the actual y coordinate.
    return y;
  }
  operator std::complex<double>() const
  {
    return std::complex<double>((x.get_d()-y.get_d()/2.),y.get_d()*M_SQRT_3_4);
  }
  void inc(int n);
  bool cont(int n);
};

extern const Eisenstein root1[6];
#endif
