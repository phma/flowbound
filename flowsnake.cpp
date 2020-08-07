/******************************************************/
/*                                                    */
/* flowsnake.cpp - flowsnake                          */
/*                                                    */
/******************************************************/
/* Copyright 2020 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include "flowsnake.h"
using namespace std;

const Eisenstein flowBase(2,-1);
const Eisenstein flowBaseConj(3,1);
vector<Segment> boundary;

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
  diff=boundary[i].b-boundary[i].a;
  lensq=diff.norm();
  cutoff=max-sqrt(4*lensq.get_d());
  for (i=0;i<boundary.size();i++)
    if (boundary[i].a.cartx()>=cutoff || boundary[i].b.cartx()>=cutoff)
      newbdy.push_back(boundary[i]);
  swap(boundary,newbdy);
}
