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
  FlowNumber a("261"),b("26.1"),c("2.61"),d(".261"),e("0.261"),f("0000000000002.61000000000000");
  FlowNumber p("14"),rp(".1111111111111111111111");
  init();
  fillTables();
  //testTables();
  for (i=0;i<42;i++)
  {
    refine();
    prune();
    cout<<i<<' '<<boundary.size()<<endl;
  }
  for (i=0;i<boundary.size();i++)
  {
    cout<<toBase7(boundary[i].a.cartx())<<','<<toBase7(boundary[i].a.carty())<<endl;
    cout<<toBase9(boundary[i].l)<<endl;
  }
  cout<<toBase7(boundary.back().b.cartx())<<','<<toBase7(boundary.back().b.carty())<<endl;
  cout<<a.toString()<<endl;
  cout<<b.toString()<<endl;
  cout<<c.toString()<<endl;
  cout<<d.toString()<<endl;
  cout<<e.toString()<<endl;
  cout<<f.toString()<<endl;
  cout<<(a+e).toString()<<endl;
  cout<<(a-e).toString()<<endl;
  cout<<(a*e).toString()<<endl;
  cout<<(p*rp).toString()<<endl;
  return 0;
}
