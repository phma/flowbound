/******************************************************/
/*                                                    */
/* flowsnake.h - flowsnake                            */
/*                                                    */
/******************************************************/
/* Copyright 2020 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include <vector>
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
