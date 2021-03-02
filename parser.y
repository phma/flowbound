/******************************************************/
/*                                                    */
/* parser.y - parser                                  */
/*                                                    */
/******************************************************/
/* Copyright 2021 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 */

%require "3.2"
%language "c++"
%define api.token.constructor
%define api.value.type variant

%{
#include "flowsnake.h"
#include "parser.h"
namespace yy
{
  parser::symbol_type yylex();
}
%}

%token
  PLUS	"+"
  MINUS	"-"
  TIMES	"*"
  OVER	"/"
  LPAR	"("
  RPAR	")"
  PREC "prec"
  INT
  FLOWNUM
;

%%

input:
  %empty
| line input
;

line:
  '\n'
| cmd '\n'
| exp '\n'
;

cmd:
  PREC INT
;

exp:
  FLOWNUM
| exp PLUS exp
;

%%

yy::parser::symbol_type yy::yylex()
{
  return 0;
}

void yy::parser::error(const std::string &msg)
{
  std::cerr<<msg<<std::endl;
}
