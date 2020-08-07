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
  for (i=0;i<33;i++)
  {
    refine();
    prune();
  }
  for (i=0;i<boundary.size();i++)
    cout<<boundary[i].a.getx()<<','<<boundary[i].a.gety()<<endl;
  return 0;
}
