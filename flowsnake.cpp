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
const Eisenstein limbBase(-25539,25807);
const FlowNumber flowOne("1"),deg60("2");
vector<Segment> boundary;
vector<uint32_t> additionTable,multiplicationTable,negationTable;
/* In these tables, the digits 0 and 1 represent themselves, but 2 is 1+ω,
 * 3 is ω, and so around the circle.
 */
const char aTable[7][7]=
{
  0,1,2,3,4,5,6,	// 0, 1, 2, 3, 4, 5, 6
  1,10,19,2,0,6,11,	// 1,13,25, 2, 0, 6,14
  2,19,18,27,3,0,1,	// 2,25,24,36, 3, 0, 1
  3,2,27,26,29,4,0,	// 3, 2,36,35,41, 4, 0
  4,0,3,29,34,37,5,	// 4, 0, 3,41,46,52, 5
  5,6,0,4,37,36,45,	// 5, 6, 0, 4,52,51,63
  6,11,1,0,5,45,44	// 6,14, 1, 0, 5,63,62
};
const char mTable[7][7]=
{
  0,0,0,0,0,0,0,
  0,1,2,3,4,5,6,
  0,2,3,4,5,6,1,
  0,3,4,5,6,1,2,
  0,4,5,6,1,2,3,
  0,5,6,1,2,3,4,
  0,6,1,2,3,4,5
};
const int pow7[]={1,7,49,343,2401,16807,117649,823543,5764801,40353607,282475249,1977326743};

void init()
{
  boundary.clear();
  boundary.push_back(Segment(Eisenstein(1,-1),Eisenstein(2,1),0));
}

void refine()
{
  vector<Segment> newbdy;
  int i;
  Eisenstein diff,p0,p1,p2,p3;
  mpz_class along;
  for (i=0;i<boundary.size();i++)
  {
    diff=(boundary[i].b-boundary[i].a)*flowBaseConj;
    p0=boundary[i].a*7;
    p3=boundary[i].b*7;
    p1=p0+diff;
    p2=p3-diff;
    along=boundary[i].l*3;
    newbdy.push_back(Segment(p0,p1,along));
    newbdy.push_back(Segment(p1,p2,along+1));
    newbdy.push_back(Segment(p2,p3,along+2));
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

string toBase9(mpz_class n)
{
  int i;
  mpz_class rem;
  string ret;
  while (n || ret.length()==0)
  {
    rem=n%9;
    n/=9;
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

uint32_t negateLimb(uint32_t a)
{
  uint32_t hi,lo;
  uint32_t ret=0;
  hi=a/pow7[6];
  lo=a%pow7[6];
  ret=negationTable[hi]*pow7[6]+negationTable[lo];
  return ret;
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
  res[0]=0;
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
      while (resDig[k]>=343 || resDig[k+1]>=343 || carry2)
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
  int i,j,neg,pos;
  for (i=0;i<343;i++)
    for (j=0;j<343;j++)
    {
      additionTable.push_back(add343(i,j));
      multiplicationTable.push_back(mul343(i,j));
    }
  additionTable.shrink_to_fit();
  multiplicationTable.shrink_to_fit();
  for (i=0;i<pow7[6];i++)
  {
    for (neg=j=0,pos=i;j<6;j++)
    {
      neg*=7;
      if (pos>=pow7[5] && pos<pow7[5]*4)
	neg+=pos/pow7[5]+3;
      if (pos>=pow7[5]*4)
	neg+=pos/pow7[5]-3;
      pos=(pos%pow7[5])*7;
    }
    negationTable.push_back(neg);
  }
  negationTable.shrink_to_fit();
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
  // Test associativity of addition
  for (i=0;i<282475249;i+=1478928)
    for (j=0;j<282475249;j+=1870697)
      for (k=0;k<282475249;k+=1578074)
      {
	intRes=addLimbs(i,j);
	left=addLimbs(intRes[0],k);
	intRes=addLimbs(j,k);
	right=addLimbs(i,intRes[0]);
	assert(left[0]==right[0]);
	assert(left[1]==right[1]);
      }
}

FlowNumber::FlowNumber()
{
  exponent=0;
}

FlowNumber::FlowNumber(std::string a)
{
  size_t i,l;
  size_t len=a.find_first_not_of("0123456.");
  if (len<a.length())
    a.erase(len);
  size_t point=a.find('.');
  if (point>a.length())
    point=a.length();
  limbs.resize((point+10)/11+(a.length()-point+9)/11);
  exponent=-((a.length()-point+9)/11);
  for (i=0;i<a.length();i++)
  {
    if (i<point)
      l=(point-i-1)/11-exponent;
    else
      l=-((i-point-1)/11)-exponent-1;
    if (a[i]!='.')
      limbs[l]=limbs[l]*7+a[i]-'0';
  }
  for (;i>point && (i-point-1)%11;i++)
    limbs[0]*=7;
  normalize();
}

string FlowNumber::toString()
{
  int i,j;
  string ret;
  for (i=limbs.size()-1;i>=0;i--)
    for (j=10;j>=0;j--)
      ret+=(char)'0'+(limbs[i]/pow7[j])%7;
  for (i=0;i<exponent;i++)
    ret+="00000000000";
  while (ret.length()/11+exponent<0)
    ret="00000000000"+ret;
  if (exponent>=0)
    ret+='.';
  else
    ret.insert(ret.length()+11*exponent,".");
  ret.erase(0,ret.find_first_not_of('0'));
  ret.erase(ret.find_last_not_of('0')+1);
  if (ret[0]=='.')
    ret="0"+ret;
  if (ret[ret.length()-1]=='.')
    ret.erase(ret.length()-1);
  return ret;
}

void FlowNumber::normalize()
{
  int i;
  for (i=0;i<limbs.size() && limbs[i]==0;i++);
  exponent+=i;
  if (i)
  {
    memmove(&limbs[0],&limbs[1],sizeof(limbs[0])*(limbs.size()-i));
    memset(&limbs[limbs.size()-i],0,sizeof(limbs[0])*i);
  }
  for (i=limbs.size()-1;i>=0 && limbs[i]==0;i--);
  limbs.resize(i+1);
  if (limbs.size()==0)
    exponent=0;
}

FlowNumber::operator complex<double>() const
{
  int i;
  Eisenstein acc;
  for (i=limbs.size()-1;i>=0;i--)
    acc=acc*limbBase+limbToEisenstein(limbs[i]);
  return (complex<double>)acc*pow((complex<double>)limbBase,exponent);
}

FlowNumber operator+(const FlowNumber &l,const FlowNumber &r)
{
  int highExp,i,j;
  FlowNumber ret;
  array<uint32_t,2> resLimb;
  ret.exponent=min(l.exponent,r.exponent);
  highExp=max(l.exponent+l.limbs.size(),r.exponent+r.limbs.size())+1;
  ret.limbs.resize(highExp-ret.exponent);
  for (i=ret.exponent;i<highExp;i++)
  {
    if (i>=l.exponent && i<(signed)(l.exponent+l.limbs.size()))
    {
      resLimb=addLimbs(ret.limbs[i-ret.exponent],l.limbs[i-l.exponent]);
      ret.limbs[i-ret.exponent]=resLimb[0];
      j=i;
      while (resLimb[1])
      {
	j++;
	resLimb=addLimbs(resLimb[1],ret.limbs[j-ret.exponent]);
	ret.limbs[j-ret.exponent]=resLimb[0];
      }
    }
    if (i>=r.exponent && i<(signed)(r.exponent+r.limbs.size()))
    {
      resLimb=addLimbs(ret.limbs[i-ret.exponent],r.limbs[i-r.exponent]);
      ret.limbs[i-ret.exponent]=resLimb[0];
      j=i;
      while (resLimb[1])
      {
	j++;
	resLimb=addLimbs(resLimb[1],ret.limbs[j-ret.exponent]);
	ret.limbs[j-ret.exponent]=resLimb[0];
      }
    }
  }
  ret.normalize();
  return ret;
}

FlowNumber operator-(const FlowNumber &l,const FlowNumber &r)
{
  int highExp,i,j;
  FlowNumber ret;
  array<uint32_t,2> resLimb;
  ret.exponent=min(l.exponent,r.exponent);
  highExp=max(l.exponent+l.limbs.size(),r.exponent+r.limbs.size())+1;
  ret.limbs.resize(highExp-ret.exponent);
  for (i=ret.exponent;i<highExp;i++)
  {
    if (i>=l.exponent && i<(signed)(l.exponent+l.limbs.size()))
    {
      resLimb=addLimbs(ret.limbs[i-ret.exponent],l.limbs[i-l.exponent]);
      ret.limbs[i-ret.exponent]=resLimb[0];
      j=i;
      while (resLimb[1])
      {
	j++;
	resLimb=addLimbs(resLimb[1],ret.limbs[j-ret.exponent]);
	ret.limbs[j-ret.exponent]=resLimb[0];
      }
    }
    if (i>=r.exponent && i<(signed)(r.exponent+r.limbs.size()))
    {
      resLimb=addLimbs(ret.limbs[i-ret.exponent],negateLimb(r.limbs[i-r.exponent]));
      ret.limbs[i-ret.exponent]=resLimb[0];
      j=i;
      while (resLimb[1])
      {
	j++;
	resLimb=addLimbs(resLimb[1],ret.limbs[j-ret.exponent]);
	ret.limbs[j-ret.exponent]=resLimb[0];
      }
    }
  }
  ret.normalize();
  return ret;
}

FlowNumber operator*(const FlowNumber &l,const FlowNumber &r)
{
  int i,j,k;
  FlowNumber ret;
  array<uint32_t,2> resLimb,prodLimb;
  ret.exponent=l.exponent+r.exponent;
  ret.limbs.resize(l.limbs.size()+r.limbs.size());
  for (i=0;i<l.limbs.size();i++)
    for (j=0;j<r.limbs.size();j++)
    {
      prodLimb=mulLimbs(l.limbs[i],r.limbs[j]);
      k=i+j;
      resLimb=addLimbs(prodLimb[0],ret.limbs[k]);
      ret.limbs[k]=resLimb[0];
      while (resLimb[1])
      {
	k++;
	resLimb=addLimbs(resLimb[1],ret.limbs[k]);
	ret.limbs[k]=resLimb[0];
      }
      k=i+j+1;
      resLimb=addLimbs(prodLimb[1],ret.limbs[k]);
      ret.limbs[k]=resLimb[0];
      while (resLimb[1])
      {
	k++;
	resLimb=addLimbs(resLimb[1],ret.limbs[k]);
	ret.limbs[k]=resLimb[0];
      }
    }
  ret.normalize();
  return ret;
}

FlowNumber operator<<(const FlowNumber &l,int n)
{
  FlowNumber ret;
  int limbShift,digShift,i;
  digShift=n%11;
  if (digShift<0)
    digShift+=11;
  limbShift=(n-digShift)/11;
  ret.limbs.resize(l.limbs.size()+1);
  for (i=0;i<l.limbs.size();i++)
  {
    ret.limbs[i]+=(l.limbs[i]%pow7[11-digShift])*pow7[digShift];
    ret.limbs[i+1]+=(l.limbs[i]/pow7[11-digShift]);
  }
  ret.exponent=l.exponent+limbShift;
  ret.normalize();
  return ret;
}

FlowNumber complexToFlowNumber(complex<double> z)
{
  int i,n,j,xp=0,sz;
  bool stop=false;
  FlowNumber inc;
  vector<FlowNumber> approx;
  vector<double> diff;
  if (z!=0.)
    xp=lrint(log(norm(z))/log(7))+2;
  n=xp;
  approx.push_back(FlowNumber());
  diff.push_back(abs((complex<double>)approx[0]-z));
  while (!stop)
  {
    inc=flowOne<<xp;
    sz=approx.size();
    for (i=0;i<sz;i++)
      for (j=0;j<6;j++)
      {
	approx.push_back(approx[i]+inc);
	diff.push_back(abs((complex<double>)approx.back()-z));
	inc=inc*deg60;
      }
    for (i=0;i<approx.size();i++)
      for (j=i-1;j>=0;j--)
	if (diff[j]>diff[j+1])
	{
	  swap(diff[j],diff[j+1]);
	  swap(approx[j],approx[j+1]);
	}
    if (diff[0]==0 || n-xp>44)
      stop=true;
    diff.resize(4);
    approx.resize(4);
    xp--;
  }
  return approx[0];
}
