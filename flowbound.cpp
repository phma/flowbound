/* Copyright 2020,2021 Pierre Abbat
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
  FlowNumber I("2.3432");
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
  cout<<a.toString()<<' '<<(complex<double>)a<<endl;
  cout<<b.toString()<<' '<<(complex<double>)b<<endl;
  cout<<c.toString()<<' '<<(complex<double>)c<<endl;
  cout<<d.toString()<<' '<<(complex<double>)d<<endl;
  cout<<e.toString()<<' '<<(complex<double>)e<<endl;
  cout<<f.toString()<<' '<<(complex<double>)f<<endl;
  cout<<(a+e).toString()<<' '<<(complex<double>)(a+e)<<endl;
  cout<<(a-e).toString()<<' '<<(complex<double>)(a-e)<<endl;
  cout<<(a*e).toString()<<' '<<(complex<double>)(a*e)<<endl;
  cout<<(p*rp).toString()<<' '<<(complex<double>)(p*rp)<<endl;
  cout<<(I).toString()<<' '<<(complex<double>)(I)<<endl;
  cout<<(I*I).toString()<<' '<<(complex<double>)(I*I)<<endl;
  return 0;
}
