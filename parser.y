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

%{
#include "flowsnake.h"
#include "parser.h"
int yylex(yy::parser::semantic_type *val);
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
