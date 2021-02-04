/******************************************************/
/*                                                    */
/* eisenstein.cpp - Eisenstein integers               */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */

#include <cstdio>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include "eisenstein.h"
//#include "ps.h"

using namespace std;
const Eisenstein root1[6] {{1,0},{1,1},{0,1},{-1,0},{-1,-1},{0,-1}};
const complex<double> omega(-0.5,M_SQRT_3_4); // this is Eisenstein(0,1)
int debugEisenstein;

mpz_class Eisenstein::numx,Eisenstein::numy,Eisenstein::denx=0,Eisenstein::deny=0,Eisenstein::quox,Eisenstein::quoy,Eisenstein::remx,Eisenstein::remy;

mpz_class _norm(mpz_class x,mpz_class y)
{
  return sqr(x)+sqr(y)-x*y;
}

Eisenstein::Eisenstein(complex<double> z)
{
  double norm0,norm1;
  y=lrint(z.imag()/M_SQRT_3_4);
  x=lrint(z.real()+y.get_d()*0.5);
  norm0=::norm(z-(complex<double> (*this)));
  y++;
  norm1=::norm(z-(complex<double> (*this)));
  if (norm1>norm0)
    y--;
  else
    norm0=norm1;
  y--,x--;
  norm1=::norm(z-(complex<double> (*this)));
  if (norm1>norm0)
    y++,x++;
  else
    norm0=norm1;
  /*x++;
  norm1=::norm(z-(complex<double> (*this)));
  if (norm1>norm0)
    x--;
  else
    norm0=norm1;*/
  y--;
  norm1=::norm(z-(complex<double> (*this)));
  if (norm1>norm0)
    y++;
  else
    norm0=norm1;
  y++,x++;
  norm1=::norm(z-(complex<double> (*this)));
  if (norm1>norm0)
    y--,x--;
  else
    norm0=norm1;
  /*x--;
  norm1=::norm(z-(complex<double> (*this)));
  if (norm1>norm0)
    x++;
  else
    norm0=norm1;*/
}

void Eisenstein::divmod(Eisenstein b)
/* Division and remainder, done together to save time
 * 1     denx       deny
 * 1+ω   denx-deny  denx
 * ω     -deny      denx-deny
 * -1    -denx      -deny
 * -1-ω  deny-denx  -denx
 * -ω    deny       deny-denx
 */
{
  int cont;
  if (this->x!=numx || this->y!=numy || b.x!=denx || b.y!=deny)
  {
    mpz_class nrm,nrm1,nd;
    numx=this->x;
    numy=this->y;
    denx=b.x;
    deny=b.y;
    nrm=b.norm();
    if (debugEisenstein)
      cout<<numx<<'+'<<numy<<"ω/"<<denx<<'+'<<deny<<"ω\n";
    // Do a rough division.
    nd=numx*denx+numy*deny-numx*deny;
    quox=round(nd.get_d()/nrm.get_d());
    nd=numy*denx-numx*deny;
    quoy=round(nd.get_d()/nrm.get_d());
    remx=numx-denx*quox+deny*quoy;
    remy=numy-denx*quoy-deny*quox+deny*quoy;
    // Adjust division so that remainder has least norm.
    // Ties are broken by < or <= for a symmetrical, but eccentric,
    // shape when dividing by LETTERMOD.
    do
    {
      cont=false; // FIXME this loop may need to be optimized
      nrm=_norm(remx,remy);
      nrm1=_norm(remx+denx-deny,remy+denx);
      if (debugEisenstein)
	cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
      if (nrm1<nrm)
      {
	remx=remx+denx-deny;
	remy=remy+denx;
	quox--;
	quoy--;
	cont-=13;
      }
      nrm=_norm(remx,remy);
      nrm1=_norm(remx+deny,remy+deny-denx);
      if (debugEisenstein)
	cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
      if (nrm1<nrm)
      {
	remx=remx+deny;
	remy=remy+deny-denx;
	quoy++;
	cont+=8;
      }
      nrm=_norm(remx,remy);
      nrm1=_norm(remx-denx,remy-deny);
      if (debugEisenstein)
	cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
      if (nrm1<nrm)
      {
	remx=remx-denx;
	remy=remy-deny;
	quox++;
	cont+=5;
      }
      nrm=_norm(remx,remy);
      nrm1=_norm(remx+deny-denx,remy-denx);
      if (debugEisenstein)
	cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
      if (nrm1<nrm)
      {
	remx=remx+deny-denx;
	remy=remy-denx;
	quox++;
	quoy++;
	cont+=13;
      }
      nrm=_norm(remx,remy);
      nrm1=_norm(remx-deny,remy+denx-deny);
      if (debugEisenstein)
	cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
      if (nrm1<nrm)
      {
	remx=remx-deny;
	remy=remy+denx-deny;
	quoy--;
	cont-=8;
      }
      nrm=_norm(remx,remy);
      nrm1=_norm(remx+denx,remy+deny);
      if (debugEisenstein)
	cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
      if (nrm1<nrm)
      {
	remx=remx+denx;
	remy=remy+deny;
	quox--;
	cont-=5;
      }
      if (debugEisenstein)
	printf("loop\n");
      } while (0);
    nrm=_norm(remx,remy);
    nrm1=_norm(remx+denx,remy+deny);
    if (debugEisenstein)
      cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
    if (nrm1<=nrm)
    {
      remx=remx+denx;
      remy=remy+deny;
      quox--;
    }
    if (debugEisenstein)
      cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω\n";
    nrm=_norm(remx,remy);
    nrm1=_norm(remx-deny+denx,remy+denx);
    if (debugEisenstein)
      cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
    if (nrm1<=nrm)
    {
      remx=remx-deny+denx;
      remy=remy+denx;
      quox--;
      quoy--;
    }
    nrm=_norm(remx,remy);
    nrm1=_norm(remx+deny,remy-denx+deny);
    if (debugEisenstein)
      cout<<"quo="<<quox<<'+'<<quoy<<"ω rem="<<remx<<'+'<<remy<<"ω nrm="<<nrm<<" nrm1="<<nrm1<<endl;
    if (nrm1<=nrm)
    {
      remx=remx+deny;
      remy=remy-denx+deny;
      quoy++;
    }
  }
}

Eisenstein Eisenstein::operator+(Eisenstein b)
{
  return Eisenstein(this->x+b.x,this->y+b.y);
}

Eisenstein& Eisenstein::operator+=(Eisenstein b)
{
  this->x+=b.x,this->y+=b.y;
  return *this;
}

Eisenstein Eisenstein::operator-()
{
  return Eisenstein(-this->x,-this->y);
}

Eisenstein Eisenstein::operator-(Eisenstein b)
{
  return Eisenstein(this->x-b.x,this->y-b.y);
}

bool operator<(const Eisenstein a,const Eisenstein b)
// These numbers are complex, so there is no consistent < operator on them.
// This operator is used only to give some order to the map.
{
  if (a.y!=b.y)
    return a.y<b.y;
  else
    return a.x<b.x;
}

Eisenstein Eisenstein::operator*(Eisenstein b)
{
  return Eisenstein(x*b.x-y*b.y,x*b.y+y*b.x-y*b.y);
}

Eisenstein& Eisenstein::operator*=(Eisenstein b)
{
  mpz_class tmp;
  tmp=x*b.x-y*b.y;
  y=x*b.y+y*b.x-y*b.y;
  x=tmp;
  return *this;
}

Eisenstein Eisenstein::operator/(Eisenstein b)
{if (b==0)
    throw(domain_error("Divide by zero Eisenstein integer"));
 divmod(b);
 return Eisenstein(quox,quoy);
 }

Eisenstein Eisenstein::operator%(Eisenstein b)
{
  if (b==0)
    return (*this); // Dividing by zero is an error, but modding by zero is not.
  else
  {
    divmod(b);
    return Eisenstein(remx,remy);
  }
}

bool Eisenstein::operator==(Eisenstein b)
{
  return this->x==b.x && this->y==b.y;
}

bool Eisenstein::operator!=(Eisenstein b)
{
  return this->x!=b.x || this->y!=b.y;
}

mpz_class Eisenstein::norm()
{
  return sqr(this->x)+sqr(this->y)-this->x*this->y;
}
