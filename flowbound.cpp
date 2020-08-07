/* Copyright 2020 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include <iostream>
#include <vector>
#include <set>
#include "eisenstein.h"
#include "flowsnake.h"
using namespace std;

int main(int argc, char *argv[])
{
  int i,j;
  init();
  for (i=0;i<42;i++)
  {
    refine();
    prune();
  }
  for (i=0;i<boundary.size();i++)
    cout<<toBase7(boundary[i].a.cartx())<<','<<toBase7(boundary[i].a.carty())<<endl;
  cout<<toBase7(boundary.back().b.cartx())<<','<<toBase7(boundary.back().b.carty())<<endl;
  return 0;
}
