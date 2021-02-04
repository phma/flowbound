/******************************************************/
/*                                                    */
/* flowsnake.cpp - flowsnake                          */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include <iostream>
#include "flowsnake.h"
using namespace std;

const Eisenstein flowBase(2,-1);
const Eisenstein flowBaseConj(3,1);
vector<Segment> boundary;
vector<uint32_t> additionTable,multiplicationTable;
/* In these tables, the digits 0 and 1 represent themselves, but 2 is 1+ω,
 * 3 is ω, and so around the circle.
 */
char aTable[7][7]=
{
  0,1,2,3,4,5,6,	// 0, 1, 2, 3, 4, 5, 6
  1,10,19,2,0,6,11,	// 1,13,25, 2, 0, 6,14
  2,19,18,27,3,0,1,	// 2,25,24,36, 3, 0, 1
  3,2,27,26,29,4,0,	// 3, 2,36,35,41, 4, 0
  4,0,3,29,34,37,5,	// 4, 0, 3,41,46,52, 5
  5,6,0,4,37,36,45,	// 5, 6, 0, 4,52,51,63
  6,11,1,0,5,45,44	// 6,14, 1, 0, 5,63,62
};
char mTable[7][7]=
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
  while (n || ret.length()==0)
  {
    rem=n%7;
    n/=7;
    ret=(char)('0'+rem.get_si())+ret;
  }
  return ret;
}

int add343(int a,int b)
{
  int aDig[3],bDig[3],resDig[4];
  int i,carry,res=0;
  for (i=0;i<3;i++)
  {
    aDig[i]=a%7;
    a/=7;
    bDig[i]=b%7;
    b/=7;
    resDig[i]=0;
  }
  resDig[3]=0;
  for (i=0;i<3;i++)
  {
    carry=0;
    resDig[i]=aTable[aDig[i]][resDig[i]];
    if (resDig[i]>6)
    {
      carry=resDig[i]/7;
      resDig[i]%=7;
    }
    resDig[i]=aTable[bDig[i]][resDig[i]];
    if (resDig[i]>6)
    { // Adding three digits cannot produce more than a two-digit number.
      carry=aTable[carry][resDig[i]/7];
      resDig[i]%=7;
    }
    resDig[i+1]=carry;
  }
  for (i=3;i>=0;i--)
    res=7*res+resDig[i];
  return res;
}

int mul343(int a,int b)
{
  int aDig[3],bDig[3],prods[3][3],resDig[6];
  int i,j,k,carry,res=0;
  for (i=0;i<3;i++)
  {
    aDig[i]=a%7;
    a/=7;
    bDig[i]=b%7;
    b/=7;
    resDig[i]=resDig[i+3]=0;
  }
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      prods[i][j]=mTable[aDig[i]][bDig[j]];
  for (i=0;i<3;i++)
  {
    for (j=0;j<3;j++)
    {
      k=i+j;
      resDig[k]=aTable[resDig[k]][prods[i][j]];
      while (resDig[k]>=7)
      {
	carry=resDig[k]/7;
	resDig[k]%=7;
	k++;
	resDig[k]=aTable[resDig[k]][carry];
      }
    }
  }
  for (i=5;i>=0;i--)
    res=7*res+resDig[i];
  return res;
}

void fillTables()
{
  int i,j;
  for (i=0;i<343;i++)
    for (j=0;j<343;j++)
    {
      additionTable.push_back(add343(i,j));
      multiplicationTable.push_back(mul343(i,j));
    }
  additionTable.shrink_to_fit();
  multiplicationTable.shrink_to_fit();
}

void testTables()
{
  int i,j;
  vector<int> realInts;
  for (i=270;i<343;i=additionTable[i*343+1])
  {
    realInts.push_back(i);
    cout<<toBase7(mpz_class(i))<<' ';
  }
  cout<<'\n'<<realInts.size()<<" integers with up to three digits\n";
}
