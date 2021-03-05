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

%token PLUS;
%token MINUS;
%token TIMES;
%token OVER;
%token LPAR;
%token RPAR;
%token PREC;
%token<VariantNumber> NUMBER;

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
  PREC NUMBER
;

exp:
  NUMBER
| exp PLUS exp
;

%%

yy::parser::symbol_type yy::yylex()
{
  int ch=std::cin.get();
  std::string buf;
  VariantNumber num;
  while (ch==' ')
    ch=std::cin.get();
  switch (ch)
  {
    case '\n':
    case EOF:
      return ch;
      break;
    case '+':
      return parser::make_PLUS();
      break;
    case '-':
      return parser::make_MINUS();
      break;
    case '*':
      return parser::make_TIMES();
      break;
    case '/':
      return parser::make_OVER();
      break;
    case '(':
      return parser::make_LPAR();
      break;
    case ')':
      return parser::make_RPAR();
      break;
    default:
      buf+=(char)ch;
      if (isalpha(ch))
      {
	while (true)
	{
	  ch=std::cin.get();
	  if (isalpha(ch))
	    buf+=(char)ch;
	  else
	    break;
	}
	std::cin.unget();
	return parser::make_PREC();
      }
      else
      {
	while (true)
	{
	  ch=std::cin.get();
	  if (isdigit(ch) || ch=='.')
	    buf+=(char)ch;
	  else
	    break;
	}
	std::cin.unget();
	num.flowNumber=FlowNumber("261");
	num.integer=8;
	num.whichValid=3;
	return parser::make_NUMBER(num);
      }
  }
}

void yy::parser::error(const std::string &msg)
{
  std::cerr<<msg<<std::endl;
}
