/******************************************************/
/*                                                    */
/* flowsnake.h - flowsnake                            */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include <vector>
#include <string>
#include "eisenstein.h"

struct Segment
{
  Eisenstein a,b;
  Segment(Eisenstein beg,Eisenstein end)
  {
    a=beg;
    b=end;
  }
  Segment()
  {
  }
  mpz_class lensq()
  {
    return (b-a).norm();
  }
};

extern std::vector<Segment> boundary;

void init();
void refine();
void prune();
std::string toBase7(mpz_class n);
void fillTables();
void testTables();

class FlowNumber
/* Complex numbers are expressed in base 2-ω with seven digits, 0 and the
 * sixth roots of 1, packed 11 to a limb.
 */
{
public:
  FlowNumber();
  FlowNumber(std::string a);
  std::string toString();
  void normalize();
private:
  std::vector<uint32_t> limbs;
  int exponent;
  friend FlowNumber operator+(const FlowNumber &l,const FlowNumber &r);
  friend FlowNumber operator-(const FlowNumber &l,const FlowNumber &r);
  friend FlowNumber operator*(const FlowNumber &l,const FlowNumber &r);
};
