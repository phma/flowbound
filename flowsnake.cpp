/******************************************************/
/*                                                    */
/* flowsnake.cpp - flowsnake                          */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include <iostream>
#include <cassert>
#include <array>
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

array<uint32_t,2> addLimbs(uint32_t a,uint32_t b)
{
  int aDig[4],bDig[4],resDig[4];
  int i,carry;
  array<uint32_t,2> res;
  for (i=0;i<4;i++)
  {
    aDig[i]=a%343;
    bDig[i]=b%343;
    a/=343;
    b/=343;
    resDig[i]=0;
  }
  for (i=0;i<4;i++)
  {
    carry=0;
    resDig[i]=additionTable[aDig[i]*343+resDig[i]];
    if (resDig[i]>342)
    {
      carry=resDig[i]/343;
      resDig[i]%=343;
    }
    resDig[i]=additionTable[bDig[i]*343+resDig[i]];
    if (resDig[i]>342)
    { // Adding three digits cannot produce more than a two-digit number.
      carry=additionTable[carry*343+resDig[i]/343];
      resDig[i]%=343;
    }
    assert(carry<343);
    resDig[i+1]=carry;
  }
  res[1]=resDig[3]/49;
  resDig[3]%=49;
  for (i=3;i>=0;i--)
    res[0]=343*res[0]+resDig[i];
  return res;
}

array<uint32_t,2> mulLimbs(uint32_t a,uint32_t b)
{
  int aDig[4],bDig[4],prods[4][4],resDig[8];
  int i,j,k,carry,carry2=0;
  array<uint32_t,2> res;
  for (i=0;i<4;i++)
  {
    aDig[i]=a%343;
    a/=343;
    bDig[i]=b%343;
    b/=343;
    resDig[i]=resDig[i+4]=0;
  }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      prods[i][j]=multiplicationTable[aDig[i]*343+bDig[j]];
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
    {
      k=i+j;
      resDig[k]=additionTable[resDig[k]*343+prods[i][j]%343];
      resDig[k+1]=additionTable[resDig[k+1]*343+prods[i][j]/343];
      while (resDig[k]>=343)
      {
	carry=additionTable[carry2*343+resDig[k]/343];
	resDig[k]%=343;
	k++;
	carry2=resDig[k]/343;
	resDig[k]=additionTable[(resDig[k]%343)*343+carry];
      }
    }
  }
  res[1]=0;
  assert(resDig[7]<7);
  for (i=7;i>=4;i--) // 76665554443 33222111000
    res[1]=343*res[1]+resDig[i];
  res[1]=7*res[1]+resDig[3]/49;
  res[0]=resDig[3]%49;
  for (i=2;i>=0;i--)
    res[0]=343*res[0]+resDig[i];
  return res;
}

Eisenstein limbToEisenstein(uint32_t n)
{
  int i,dig;
  Eisenstein ret,powBase=1;
  while (n)
  {
    dig=n%7;
    n/=7;
    if (dig)
      ret+=powBase*root1[dig-1];
    powBase*=flowBase;
  }
  return ret;
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
  int i,j,k;
  vector<int> realInts,upOneInts;
  array<uint32_t,2> intRes,left,right;
  for (i=270;i<343;i=additionTable[i*343+1])
  {
    realInts.push_back(i);
    cout<<toBase7(mpz_class(i))<<' ';
  }
  cout<<'\n';
  for (i=269;i<343;i=additionTable[i*343+1])
  {
    upOneInts.push_back(i);
    cout<<toBase7(mpz_class(i))<<' ';
  }
  cout<<'\n'<<realInts.size()<<" integers with up to three digits\n";
  for (i=0;i<realInts.size();i++)
    for (j=0;j<upOneInts.size();j++)
      assert(limbToEisenstein(multiplicationTable[realInts[i]*343+upOneInts[j]])==
	     limbToEisenstein(realInts[i])*limbToEisenstein(upOneInts[j]));
  // Test associativity of multiplication with intermediate results up to 7**11
  assert(mulLimbs(2667,2736)[0]==mulLimbs(2736,2667)[0]);
  assert(mulLimbs(16260,2736)[0]==mulLimbs(2736,16260)[0]);
  assert(mulLimbs(2667,39302898)[0]==mulLimbs(39302898,2667)[0]);
  assert(mulLimbs(16260,11137910)[0]==mulLimbs(11137910,16260)[0]);
  mulLimbs(5472,15718);
  for (i=0;i<16384;i+=381)
    for (j=0;j<117649;j+=2736)
      for (k=0;k<16807;k+=542)
      {
	intRes=mulLimbs(i,j);
	left=mulLimbs(intRes[0],k);
	intRes=mulLimbs(j,k);
	right=mulLimbs(i,intRes[0]);
	assert(left[0]==right[0]);
	assert(left[1]==right[1]);
      }
}
