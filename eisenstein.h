/******************************************************/
/*                                                    */
/* eisenstein.h - Eisenstein integers                 */
/*                                                    */
/******************************************************/
/* Copyright 2020 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#ifndef EISENSTEIN_H
#define EISENSTEIN_H

#include <cmath>
#include <cstdlib>
#include <map>
#include <complex>
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
  int x,y; // x is the real part, y is at 120°
  static int numx,numy,denx,deny,quox,quoy,remx,remy;
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
  unsigned long norm();
  int pageinx(int size,int nelts);
  int pageinx();
  int letterinx();
  int getx()
  {
    return x;
  }
  int gety()
  {
    return y;
  }
  operator std::complex<double>() const
  {
    return std::complex<double>(x-y/2.,y*M_SQRT_3_4);
  }
  void inc(int n);
  bool cont(int n);
};

Eisenstein start(int n);
Eisenstein nthEisenstein(int n,int size,int nelts);
extern const Eisenstein LETTERMOD,PAGEMOD;
extern int debugEisenstein;

template <typename T> class harray
{
  std::map<Eisenstein,T *> index;
public:
  T& operator[](Eisenstein i);
  void clear();
};

template <typename T> T& harray<T>::operator[](Eisenstein i)
{
  Eisenstein q,r;
  q=i/PAGEMOD;
  r=i%PAGEMOD;
  if (!index[q])
    index[q]=(T*)calloc(PAGESIZE,sizeof(T));
  return index[q][r.pageinx()];
}

template <typename T> void harray<T>::clear()
{
  typename std::map<Eisenstein,T *>::iterator i;
  for (i=index.start();i!=index.end();i++)
  {
    free(i->second);
    i->second=NULL;
  }
}

int region(std::complex<double> z);

extern harray<char> hletters,hbits;
#endif
