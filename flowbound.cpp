/* Copyright 2020,2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */
#include <iostream>
#include <vector>
#include <set>
#include "eisenstein.h"
#include "flowsnake.h"
#include "parser.h"
using namespace std;

int main(int argc, char *argv[])
{
  int i,j;
  FlowNumber a("261"),b("26.1"),c("2.61"),d(".261"),e("0.261"),f("0000000000002.61000000000000");
  FlowNumber g("200000000000"),h(".000000000001");
  FlowNumber p("14"),rp(".1111111111111111111111");
  FlowNumber imag,half,frac,east;
  FlowNumber one("1");
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
  imag=complexToFlowNumber(complex<double>(0,1));
  half=complexToFlowNumber(complex<double>(0.5,0));
  frac=complexToFlowNumber(complex<double>(1/8.,0));
  east=complexToFlowNumber(complex<double>(.5749186263504634,.07052452345448619));
  /* east is 1/3 what's output by the above, because init sets the segment ends
   * to Eisenstein integers which are 3 times the corners of the boundary.
   */
  cout<<a.toString()<<' '<<(complex<double>)a<<endl;
  cout<<b.toString()<<' '<<(complex<double>)b<<endl;
  cout<<c.toString()<<' '<<(complex<double>)c<<endl;
  cout<<d.toString()<<' '<<(complex<double>)d<<endl;
  cout<<e.toString()<<' '<<(complex<double>)e<<endl;
  cout<<f.toString()<<' '<<(complex<double>)f<<endl;
  cout<<(a+e).toString()<<' '<<(complex<double>)(a+e)<<endl;
  cout<<(g+h).toString()<<' '<<(complex<double>)(g+h)<<endl;
  cout<<(a-e).toString()<<' '<<(complex<double>)(a-e)<<endl;
  cout<<(a*e).toString()<<' '<<(complex<double>)(a*e)<<endl;
  cout<<(p*rp).toString()<<' '<<(complex<double>)(p*rp)<<endl;
  cout<<(a/a).toString()<<' '<<(complex<double>)(a/a)<<endl;
  cout<<(a/d).toString()<<' '<<(complex<double>)(a/d)<<endl;
  cout<<(d/a).toString()<<' '<<(complex<double>)(d/a)<<endl;
  cout<<(one/imag).toString()<<' '<<(complex<double>)(one/imag)<<endl;
  cout<<(imag).toString()<<' '<<(complex<double>)(imag)<<endl;
  cout<<(imag*imag).toString()<<' '<<(complex<double>)(imag*imag)<<endl;
  cout<<(half).toString()<<' '<<(complex<double>)(half)<<endl;
  cout<<(frac).toString()<<' '<<(complex<double>)(frac)<<endl;
  cout<<(east).toString()<<' '<<(complex<double>)(east)<<endl;
  return 0;
}
