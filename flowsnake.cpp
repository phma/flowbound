/******************************************************/
/*                                                    */
/* flowsnake.cpp - flowsnake                          */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include "flowsnake.h"
using namespace std;

const Eisenstein flowBase(2,-1);
const Eisenstein flowBaseConj(3,1);
vector<Segment> boundary;
vector<uint32_t> additionTable,multiplicationTable;
/* In these tables, the digits 0 and 1 represent themselves, but 2 is 1+ω,
 * 3 is ω, and so around the circle.
 */
char aTable[49]=
{
  0,1,2,3,4,5,6,	// 0, 1, 2, 3, 4, 5, 6
  1,10,19,2,0,6,11,	// 1,13,25, 2, 0, 6,14
  2,19,18,27,3,0,1,	// 2,25,24,36, 3, 0, 1
  3,2,27,26,29,4,0,	// 3, 2,36,35,41, 4, 0
  4,0,3,29,34,37,5,	// 4, 0, 3,41,46,52, 5
  5,6,0,4,37,36,45,	// 5, 6, 0, 4,52,51,63
  6,11,1,0,5,45,44	// 6,14, 1, 0, 5,63,62
};
char mTable[49]=
{
  0,0,0,0,0,0,0,
  0,1,2,3,4,5,6,
  0,2,3,4,5,6,1,
  0,3,4,5,6,1,2,
  0,4,5,6,1,2,3,
  0,5,6,1,2,3,4,
  0,6,1,2,3,4,5
};

void init()
{
  boundary.clear();
  boundary.push_back(Segment(Eisenstein(1,-1),Eisenstein(2,1)));
}

void refine()
{
  vector<Segment> newbdy;
  int i;
  Eisenstein diff,p0,p1,p2,p3;
  for (i=0;i<boundary.size();i++)
  {
    diff=(boundary[i].b-boundary[i].a)*flowBaseConj;
    p0=boundary[i].a*7;
    p3=boundary[i].b*7;
    p1=p0+diff;
    p2=p3-diff;
    newbdy.push_back(Segment(p0,p1));
    newbdy.push_back(Segment(p1,p2));
    newbdy.push_back(Segment(p2,p3));
  }
  swap(boundary,newbdy);
}

void prune()
{
  vector<Segment> newbdy;
  int i;
  Eisenstein diff;
  mpz_class max,lensq,cutoff;
  for (i=0;i<boundary.size();i++)
  {
    if (boundary[i].a.cartx()>max)
      max=boundary[i].a.cartx();
    if (boundary[i].a.cartx()>max)
      max=boundary[i].a.cartx();
  }
  diff=boundary[0].b-boundary[0].a;
  lensq=diff.norm();
  cutoff=max-sqrt(4*lensq.get_d());
  for (i=0;i<boundary.size();i++)
    if (boundary[i].a.cartx()>=cutoff || boundary[i].b.cartx()>=cutoff)
      newbdy.push_back(boundary[i]);
  swap(boundary,newbdy);
}

string toBase7(mpz_class n)
{
  int i;
  mpz_class rem;
  string ret;
  while (n)
  {
    rem=n%7;
    n/=7;
    ret=(char)('0'+rem.get_si())+ret;
  }
  return ret;
}
