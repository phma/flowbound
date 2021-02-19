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
#include <array>
#include "eisenstein.h"

extern const Eisenstein flowBase;

struct Segment
{
  Eisenstein a,b;
  mpz_class l;
  Segment(Eisenstein beg,Eisenstein end,mpz_class along)
  {
    a=beg;
    b=end;
    l=along;
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
std::string toBase9(mpz_class n);
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
  void setPrecision(int p,bool t); // t is true for relative
  std::string toString();
  void normalize();
  operator std::complex<double>() const;
  std::array<int,3> msd() const;
private:
  static int precision;
  static bool precType;
  std::vector<uint32_t> limbs;
  int exponent;
  friend FlowNumber operator+(const FlowNumber &l,const FlowNumber &r);
  friend FlowNumber operator-(const FlowNumber &l,const FlowNumber &r);
  friend FlowNumber operator*(const FlowNumber &l,const FlowNumber &r);
  friend FlowNumber operator/(FlowNumber l,const FlowNumber &r);
  friend FlowNumber operator<<(const FlowNumber &l,int n);
};

FlowNumber complexToFlowNumber(std::complex<double> z);
