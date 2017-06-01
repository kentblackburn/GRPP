
/*  A Bison parser, made from grpp_yacc.y
    by GNU Bison version 1.28  */

#define YYBISON 1  /* Identify Bison output.  */

#define	NL	257
#define	SEMI	258
#define	EQ	259
#define	SYMEQ	260
#define	MARKER	261
#define	NUMBER	262
#define	VARIABLE	263
#define	FUNCTION	264
#define	TYPE	265
#define	RANK1	266
#define	RANK2	267
#define	RANK3	268
#define	RANK4	269
#define	COORD	270
#define	INDICE	271
#define	ARRAY	272
#define	VECTOR	273
#define	ANTISYM	274
#define	COMMENT	275
#define	COVDERIV	276
#define	LIEDERIV	277
#define	UMINUS	278

#line 1 "grpp_yacc.y"


/*
	File:		grpp_yacc.y
	Description:	Yacc Parser for GRPP
	Author:		James Kent Blackburn
	Date:		Jan 2000

	Copyright (c) 1992-2000 
	By James Kent Blackburn
	All Rights Reserved
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define		TRUE	1
#define		FALSE	0

char		msg[80] = "";

extern int	coord_flag, index_flag;
extern int	n_dim,n_ind,no_lines_in,no_lines_out,gr_debug,line_flg;

extern char     *infile,*hdrfile;
extern char	coord_list[50],
		upper_coord[25],
		lower_coord[25],
		index_list[50],
		upper_index[25],
		lower_index[25];

int	is_vector(char*);
int	is_index(char*);
int	is_coord(char*);
void	write_line_no(void);
void	tensor_structures(void);
void	tensor_declare(int, char*, char*);
void	total_antisymmetric_tensor(char*);
void	covariant_derivative(char*, char*, char*);
void	lie_derivative(char*, char*, char*);
int	permutation_sort(int, int*);
void	vector_equation(char*, char*, char*);
void	equation_vector(char*, char*, char*);
void	tensor_eq_array(int, char*, char*, char*);
void	array_eq_expression(char*, char*, char*);
void	array_elements(char*, char*, int, int, int, int, int);
void	tensor_eq_expression(int, char*, char*, char*);
void	symten_eq_expression(int, char*, char*, char*);
void	expression_replace(char, char, char*);
void	remove_dollars(char*);
void	negate_expression(char *);
void	print_equation(char*, char*, char*);
void	expand_dummy_indices(int, char*);
void	yyerror(char*);
int	sum_expressions(char*, char*, char, char*);
int	tensor_product(char*, char*, char*);
void	get_tensor_name(char*,char*);
int	get_rank(char*);
int	get_all_indices(char*, char*);
int	get_free_indices(char*, char*);
int	get_dummy_indices(char*, char*);
int	cmp_all_indices(char*, char*);
int	cmp_free_indices(char*, char*);
int	parse_tensor(int, char*, char*, char*);
int	chcomp(char*, char*);
extern void *qsort(void *base, size_t nmemb, size_t size,
		   int (*compar)(char *, char *));


#line 73 "grpp_yacc.y"
typedef union	{
	 char	letter;
	 char	string[8192];
	} YYSTYPE;
#include <stdio.h>

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif



#define	YYFINAL		184
#define	YYFLAG		-32768
#define	YYNTBASE	33

#define YYTRANSLATE(x) ((unsigned)(x) <= 278 ? yytranslate[x] : 38)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,    30,
    31,    26,    24,    29,    25,     2,    27,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    32,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     3,     4,     5,     6,
     7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    28
};

#if YYDEBUG != 0
static const short yyprhs[] = {     0,
     0,     1,     4,     6,     8,    10,    12,    14,    16,    18,
    22,    26,    30,    34,    40,    46,    52,    58,    66,    74,
    82,    90,   100,   110,   120,   130,   142,   154,   166,   178,
   192,   206,   220,   234,   240,   243,   248,   253,   258,   263,
   268,   273,   278,   283,   288,   293,   298,   303,   308,   313,
   318,   320,   322,   324,   326,   328,   330,   332,   339,   346,
   353,   360,   367,   374,   381,   388,   392,   396,   400,   404,
   407
};

static const short yyrhs[] = {    -1,
    33,    34,     0,     3,     0,     7,     0,    21,     0,    16,
     0,    17,     0,    35,     0,    36,     0,    11,    12,     4,
     0,    11,    13,     4,     0,    11,    14,     4,     0,    11,
    15,     4,     0,    11,    12,    29,    12,     4,     0,    11,
    13,    29,    13,     4,     0,    11,    14,    29,    14,     4,
     0,    11,    15,    29,    15,     4,     0,    11,    12,    29,
    12,    29,    12,     4,     0,    11,    13,    29,    13,    29,
    13,     4,     0,    11,    14,    29,    14,    29,    14,     4,
     0,    11,    15,    29,    15,    29,    15,     4,     0,    11,
    12,    29,    12,    29,    12,    29,    12,     4,     0,    11,
    13,    29,    13,    29,    13,    29,    13,     4,     0,    11,
    14,    29,    14,    29,    14,    29,    14,     4,     0,    11,
    15,    29,    15,    29,    15,    29,    15,     4,     0,    11,
    12,    29,    12,    29,    12,    29,    12,    29,    12,     4,
     0,    11,    13,    29,    13,    29,    13,    29,    13,    29,
    13,     4,     0,    11,    14,    29,    14,    29,    14,    29,
    14,    29,    14,     4,     0,    11,    15,    29,    15,    29,
    15,    29,    15,    29,    15,     4,     0,    11,    12,    29,
    12,    29,    12,    29,    12,    29,    12,    29,    12,     4,
     0,    11,    13,    29,    13,    29,    13,    29,    13,    29,
    13,    29,    13,     4,     0,    11,    14,    29,    14,    29,
    14,    29,    14,    29,    14,    29,    14,     4,     0,    11,
    15,    29,    15,    29,    15,    29,    15,    29,    15,    29,
    15,     4,     0,    20,    30,     9,    31,     4,     0,    10,
     4,     0,    12,     5,    19,     4,     0,    19,     5,    12,
     4,     0,    12,     5,    18,     4,     0,    13,     5,    18,
     4,     0,    14,     5,    18,     4,     0,    15,     5,    18,
     4,     0,    18,     5,    37,     4,     0,     9,     5,    37,
     4,     0,    12,     5,    37,     4,     0,    13,     5,    37,
     4,     0,    14,     5,    37,     4,     0,    15,     5,    37,
     4,     0,     9,     6,    37,     4,     0,    12,     6,    37,
     4,     0,    13,     6,    37,     4,     0,     8,     0,     9,
     0,    10,     0,    12,     0,    13,     0,    14,     0,    15,
     0,    22,    30,    12,    32,    14,    31,     0,    22,    30,
    13,    32,    14,    31,     0,    22,    30,    14,    32,    14,
    31,     0,    22,    30,    15,    32,    14,    31,     0,    23,
    30,    12,    32,    12,    31,     0,    23,    30,    12,    32,
    13,    31,     0,    23,    30,    12,    32,    14,    31,     0,
    23,    30,    12,    32,    15,    31,     0,    37,    24,    37,
     0,    37,    25,    37,     0,    37,    26,    37,     0,    37,
    27,    37,     0,    25,    37,     0,    30,    37,    31,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
    92,    93,    96,    99,   101,   103,   109,   114,   115,   118,
   120,   122,   124,   126,   129,   132,   135,   138,   142,   146,
   150,   154,   159,   164,   169,   174,   180,   186,   192,   198,
   205,   212,   219,   228,   230,   232,   234,   236,   238,   240,
   242,   244,   246,   248,   250,   252,   254,   256,   258,   260,
   264,   268,   272,   276,   281,   286,   291,   296,   298,   300,
   302,   304,   306,   308,   310,   312,   323,   334,   345,   351,
   356
};
#endif


#if YYDEBUG != 0 || defined (YYERROR_VERBOSE)

static const char * const yytname[] = {   "$","error","$undefined.","NL","SEMI",
"EQ","SYMEQ","MARKER","NUMBER","VARIABLE","FUNCTION","TYPE","RANK1","RANK2",
"RANK3","RANK4","COORD","INDICE","ARRAY","VECTOR","ANTISYM","COMMENT","COVDERIV",
"LIEDERIV","'+'","'-'","'*'","'/'","UMINUS","','","'('","')'","':'","lines",
"line","declare","assign","expression", NULL
};
#endif

static const short yyr1[] = {     0,
    33,    33,    34,    34,    34,    34,    34,    34,    34,    35,
    35,    35,    35,    35,    35,    35,    35,    35,    35,    35,
    35,    35,    35,    35,    35,    35,    35,    35,    35,    35,
    35,    35,    35,    36,    36,    36,    36,    36,    36,    36,
    36,    36,    36,    36,    36,    36,    36,    36,    36,    36,
    37,    37,    37,    37,    37,    37,    37,    37,    37,    37,
    37,    37,    37,    37,    37,    37,    37,    37,    37,    37,
    37
};

static const short yyr2[] = {     0,
     0,     2,     1,     1,     1,     1,     1,     1,     1,     3,
     3,     3,     3,     5,     5,     5,     5,     7,     7,     7,
     7,     9,     9,     9,     9,    11,    11,    11,    11,    13,
    13,    13,    13,     5,     2,     4,     4,     4,     4,     4,
     4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
     1,     1,     1,     1,     1,     1,     1,     6,     6,     6,
     6,     6,     6,     6,     6,     3,     3,     3,     3,     2,
     3
};

static const short yydefact[] = {     1,
     0,     3,     4,     0,     0,     0,     0,     0,     0,     0,
     6,     7,     0,     0,     0,     5,     2,     8,     9,     0,
     0,    35,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,    51,    52,    53,    54,    55,
    56,    57,     0,     0,     0,     0,     0,     0,    10,     0,
    11,     0,    12,     0,    13,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,    70,     0,    43,     0,     0,     0,     0,    48,
     0,     0,     0,     0,    38,    36,    44,    49,    39,    45,
    50,    40,    46,    41,    47,    42,    37,     0,     0,     0,
     0,     0,     0,    71,    66,    67,    68,    69,    14,     0,
    15,     0,    16,     0,    17,     0,    34,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,    18,     0,    19,     0,    20,     0,
    21,     0,    58,    59,    60,    61,    62,    63,    64,    65,
     0,     0,     0,     0,    22,     0,    23,     0,    24,     0,
    25,     0,     0,     0,     0,     0,    26,     0,    27,     0,
    28,     0,    29,     0,     0,     0,     0,     0,    30,    31,
    32,    33,     0,     0
};

static const short yydefgoto[] = {     1,
    17,    18,    19,    47
};

static const short yypact[] = {-32768,
   221,-32768,-32768,    60,    46,    94,    70,    74,   105,   114,
-32768,-32768,   115,   116,    34,-32768,-32768,-32768,-32768,   197,
   197,-32768,    -3,    -2,    -1,     0,   121,   197,   140,   197,
   159,   178,   197,    56,    65,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,    48,    52,   197,   197,    27,    45,-32768,    80,
-32768,   110,-32768,   118,-32768,   109,   133,   134,    59,    63,
   137,    69,    73,   138,    77,   141,    87,   101,   143,   125,
   201,   145,-32768,    91,-32768,   197,   197,   197,   197,-32768,
     1,     9,    10,    11,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,   155,   128,   129,
   132,   144,   146,-32768,     6,     6,     6,     6,-32768,   154,
-32768,   162,-32768,   165,-32768,   168,-32768,   166,   171,   180,
   181,   231,    12,    13,    14,    15,   167,   173,   186,   187,
   192,   194,   195,   198,-32768,   185,-32768,   189,-32768,   233,
-32768,   184,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
    16,    17,    18,    19,-32768,   236,-32768,   237,-32768,   235,
-32768,   238,    30,    31,    32,    33,-32768,   239,-32768,   241,
-32768,   242,-32768,   240,   248,   253,   254,   255,-32768,-32768,
-32768,-32768,   260,-32768
};

static const short yypgoto[] = {-32768,
-32768,-32768,-32768,   -21
};


#define	YYLAST		260


static const short yytable[] = {    48,
    49,    51,    53,    55,   109,    59,    60,    62,    63,    65,
    67,    68,   111,   113,   115,   135,   137,   139,   141,   155,
   157,   159,   161,    73,    74,    50,    52,    54,    56,   110,
    75,    78,    79,   167,   169,   171,   173,   112,   114,   116,
   136,   138,   140,   142,   156,   158,   160,   162,    80,    22,
    76,    77,    78,    79,   105,   106,   107,   108,   168,   170,
   172,   174,    87,    35,    20,    21,    88,    69,    76,    77,
    78,    79,    90,    70,    27,    28,    91,    71,    29,    30,
    93,    72,    76,    77,    78,    79,    76,    77,    78,    79,
    95,    81,    76,    77,    78,    79,    76,    77,    78,    79,
    76,    77,    78,    79,    96,    23,    24,    25,    26,    31,
    76,    77,    78,    79,    76,    77,    78,    79,    32,    33,
    34,   104,    82,    84,    76,    77,    78,    79,    36,    37,
    38,    83,    39,    40,    41,    42,    85,    86,    57,    58,
    89,    92,    43,    44,    94,    45,    97,    36,    37,    38,
    46,    39,    40,    41,    42,    98,   103,    61,   117,   118,
   119,    43,    44,   120,    45,   123,    36,    37,    38,    46,
    39,    40,    41,    42,   124,   121,    64,   122,   125,   127,
    43,    44,   126,    45,   128,    36,    37,    38,    46,    39,
    40,    41,    42,   129,   130,    66,   151,   143,   154,    43,
    44,   152,    45,   144,    36,    37,    38,    46,    39,    40,
    41,    42,    99,   100,   101,   102,   145,   146,    43,    44,
   183,    45,   147,     2,   148,   149,    46,     3,   150,     4,
     5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
    15,    16,   131,   132,   133,   134,   153,   163,   165,   164,
   175,   179,   166,   176,   178,   177,   180,   181,   182,   184
};

static const short yycheck[] = {    21,
     4,     4,     4,     4,     4,    27,    28,    29,    30,    31,
    32,    33,     4,     4,     4,     4,     4,     4,     4,     4,
     4,     4,     4,    45,    46,    29,    29,    29,    29,    29,
     4,    26,    27,     4,     4,     4,     4,    29,    29,    29,
    29,    29,    29,    29,    29,    29,    29,    29,     4,     4,
    24,    25,    26,    27,    76,    77,    78,    79,    29,    29,
    29,    29,     4,    30,     5,     6,     4,    12,    24,    25,
    26,    27,     4,     9,     5,     6,     4,    30,     5,     6,
     4,    30,    24,    25,    26,    27,    24,    25,    26,    27,
     4,    12,    24,    25,    26,    27,    24,    25,    26,    27,
    24,    25,    26,    27,     4,    12,    13,    14,    15,     5,
    24,    25,    26,    27,    24,    25,    26,    27,     5,     5,
     5,    31,    13,    15,    24,    25,    26,    27,     8,     9,
    10,    14,    12,    13,    14,    15,     4,     4,    18,    19,
     4,     4,    22,    23,     4,    25,     4,     8,     9,    10,
    30,    12,    13,    14,    15,    31,    12,    18,     4,    32,
    32,    22,    23,    32,    25,    12,     8,     9,    10,    30,
    12,    13,    14,    15,    13,    32,    18,    32,    14,    14,
    22,    23,    15,    25,    14,     8,     9,    10,    30,    12,
    13,    14,    15,    14,    14,    18,    12,    31,    15,    22,
    23,    13,    25,    31,     8,     9,    10,    30,    12,    13,
    14,    15,    12,    13,    14,    15,    31,    31,    22,    23,
     0,    25,    31,     3,    31,    31,    30,     7,    31,     9,
    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
    20,    21,    12,    13,    14,    15,    14,    12,    14,    13,
    12,     4,    15,    13,    15,    14,     4,     4,     4,     0
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/lib/bison.simple"
/* This file comes from bison-1.28.  */

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

#ifndef YYSTACK_USE_ALLOCA
#ifdef alloca
#define YYSTACK_USE_ALLOCA
#else /* alloca not defined */
#ifdef __GNUC__
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi) || (defined (__sun) && defined (__i386))
#define YYSTACK_USE_ALLOCA
#include <alloca.h>
#else /* not sparc */
/* We think this test detects Watcom and Microsoft C.  */
/* This used to test MSDOS, but that is a bad idea
   since that symbol is in the user namespace.  */
#if (defined (_MSDOS) || defined (_MSDOS_)) && !defined (__TURBOC__)
#if 0 /* No need for malloc.h, which pollutes the namespace;
	 instead, just don't use alloca.  */
#include <malloc.h>
#endif
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
/* I don't know what this was needed for, but it pollutes the namespace.
   So I turned it off.   rms, 2 May 1997.  */
/* #include <malloc.h>  */
 #pragma alloca
#define YYSTACK_USE_ALLOCA
#else /* not MSDOS, or __TURBOC__, or _AIX */
#if 0
#ifdef __hpux /* haible@ilog.fr says this works for HPUX 9.05 and up,
		 and on HPUX 10.  Eventually we can turn this on.  */
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#endif /* __hpux */
#endif
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not sparc */
#endif /* not GNU C */
#endif /* alloca not defined */
#endif /* YYSTACK_USE_ALLOCA not defined */

#ifdef YYSTACK_USE_ALLOCA
#define YYSTACK_ALLOC alloca
#else
#define YYSTACK_ALLOC malloc
#endif

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    { yychar = (token), yylval = (value);			\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { yyerror ("syntax error: cannot back up"); YYERROR; }	\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

#ifndef YYPURE
#define YYLEX		yylex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, &yylloc, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval, &yylloc)
#endif
#else /* not YYLSP_NEEDED */
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval)
#endif
#endif /* not YYLSP_NEEDED */
#endif

/* If nonreentrant, generate the variables here */

#ifndef YYPURE

int	yychar;			/*  the lookahead symbol		*/
YYSTYPE	yylval;			/*  the semantic value of the		*/
				/*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE yylloc;			/*  location data for the lookahead	*/
				/*  symbol				*/
#endif

int yynerrs;			/*  number of parse errors so far       */
#endif  /* not YYPURE */

#if YYDEBUG != 0
int yydebug;			/*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef	YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

/* Define __yy_memcpy.  Note that the size argument
   should be passed with type unsigned int, because that is what the non-GCC
   definitions require.  With GCC, __builtin_memcpy takes an arg
   of type size_t, but it can handle unsigned int.  */

#if __GNUC__ > 1		/* GNU C and GNU C++ define this.  */
#define __yy_memcpy(TO,FROM,COUNT)	__builtin_memcpy(TO,FROM,COUNT)
#else				/* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (to, from, count)
     char *to;
     char *from;
     unsigned int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (char *to, char *from, unsigned int count)
{
  register char *t = to;
  register char *f = from;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

#line 217 "/usr/lib/bison.simple"

/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
#ifdef __cplusplus
#define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else /* not __cplusplus */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
#endif /* not __cplusplus */
#else /* not YYPARSE_PARAM */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif /* not YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
#ifdef YYPARSE_PARAM
int yyparse (void *);
#else
int yyparse (void);
#endif
#endif

int
yyparse(YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YYSTYPE *yyvsp;
  int yyerrstatus;	/*  number of tokens to shift before error messages enabled */
  int yychar1 = 0;		/*  lookahead token as an internal (translated) token number */

  short	yyssa[YYINITDEPTH];	/*  the state stack			*/
  YYSTYPE yyvsa[YYINITDEPTH];	/*  the semantic value stack		*/

  short *yyss = yyssa;		/*  refer to the stacks thru separate pointers */
  YYSTYPE *yyvs = yyvsa;	/*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE yylsa[YYINITDEPTH];	/*  the location stack			*/
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;
  int yyfree_stacks = 0;

#ifdef YYPURE
  int yychar;
  YYSTYPE yylval;
  int yynerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE yylloc;
#endif
#endif

  YYSTYPE yyval;		/*  the variable used to return		*/
				/*  semantic values from the action	*/
				/*  routines				*/

  int yylen;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Starting parse\n");
#endif

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YYLSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
yynewstate:

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YYLSP_NEEDED
      YYLTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YYLSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YYLSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  yyerror("parser stack overflow");
	  if (yyfree_stacks)
	    {
	      free (yyss);
	      free (yyvs);
#ifdef YYLSP_NEEDED
	      free (yyls);
#endif
	    }
	  return 2;
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
#ifndef YYSTACK_USE_ALLOCA
      yyfree_stacks = 1;
#endif
      yyss = (short *) YYSTACK_ALLOC (yystacksize * sizeof (*yyssp));
      __yy_memcpy ((char *)yyss, (char *)yyss1,
		   size * (unsigned int) sizeof (*yyssp));
      yyvs = (YYSTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yyvsp));
      __yy_memcpy ((char *)yyvs, (char *)yyvs1,
		   size * (unsigned int) sizeof (*yyvsp));
#ifdef YYLSP_NEEDED
      yyls = (YYLTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yylsp));
      __yy_memcpy ((char *)yyls, (char *)yyls1,
		   size * (unsigned int) sizeof (*yylsp));
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YYLSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  goto yybackup;
 yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Reading a token: ");
#endif
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(yychar);

#if YYDEBUG != 0
      if (yydebug)
	{
	  fprintf (stderr, "Next token is %d (%s", yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting token %d (%s), ", yychar, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  goto yynewstate;

/* Do the default action for the current state.  */
yydefault:

  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
yyreduce:
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (yydebug)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


  switch (yyn) {

case 3:
#line 97 "grpp_yacc.y"
{no_lines_out++; printf("\n");
			 write_line_no();;
    break;}
case 4:
#line 100 "grpp_yacc.y"
{printf("%s",yyvsp[0].string);;
    break;}
case 5:
#line 102 "grpp_yacc.y"
{printf("%s",yyvsp[0].string);;
    break;}
case 6:
#line 104 "grpp_yacc.y"
{printf("/* COORDINATES: %s */\n",coord_list);
			 printf("/* DIMENSION: %d */\n",n_dim);
			 printf("\n#include \"%s\" \n",hdrfile);
			 no_lines_out += 4;
			 tensor_structures();;
    break;}
case 7:
#line 110 "grpp_yacc.y"
{printf("/* INDICES: %s */\n",index_list);
			 printf("/* TOTAL: %d */\n",n_ind);
			 no_lines_out += 2;
			 tensor_structures();;
    break;}
case 10:
#line 119 "grpp_yacc.y"
{tensor_declare(1,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 11:
#line 121 "grpp_yacc.y"
{tensor_declare(2,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 12:
#line 123 "grpp_yacc.y"
{tensor_declare(3,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 13:
#line 125 "grpp_yacc.y"
{tensor_declare(4,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 14:
#line 127 "grpp_yacc.y"
{tensor_declare(1,yyvsp[-4].string,yyvsp[-3].string);
			 tensor_declare(1,yyvsp[-4].string,yyvsp[-1].string);;
    break;}
case 15:
#line 130 "grpp_yacc.y"
{tensor_declare(2,yyvsp[-4].string,yyvsp[-3].string);
			 tensor_declare(2,yyvsp[-4].string,yyvsp[-1].string);;
    break;}
case 16:
#line 133 "grpp_yacc.y"
{tensor_declare(3,yyvsp[-4].string,yyvsp[-3].string);
			 tensor_declare(3,yyvsp[-4].string,yyvsp[-1].string);;
    break;}
case 17:
#line 136 "grpp_yacc.y"
{tensor_declare(4,yyvsp[-4].string,yyvsp[-3].string);
			 tensor_declare(4,yyvsp[-4].string,yyvsp[-1].string);;
    break;}
case 18:
#line 139 "grpp_yacc.y"
{tensor_declare(1,yyvsp[-6].string,yyvsp[-5].string);
			 tensor_declare(1,yyvsp[-6].string,yyvsp[-3].string);
			 tensor_declare(1,yyvsp[-6].string,yyvsp[-1].string);;
    break;}
case 19:
#line 143 "grpp_yacc.y"
{tensor_declare(2,yyvsp[-6].string,yyvsp[-5].string);
			 tensor_declare(2,yyvsp[-6].string,yyvsp[-3].string);
			 tensor_declare(2,yyvsp[-6].string,yyvsp[-1].string);;
    break;}
case 20:
#line 147 "grpp_yacc.y"
{tensor_declare(3,yyvsp[-6].string,yyvsp[-5].string);
			 tensor_declare(3,yyvsp[-6].string,yyvsp[-3].string);
			 tensor_declare(3,yyvsp[-6].string,yyvsp[-1].string);;
    break;}
case 21:
#line 151 "grpp_yacc.y"
{tensor_declare(4,yyvsp[-6].string,yyvsp[-5].string);
			 tensor_declare(4,yyvsp[-6].string,yyvsp[-3].string);
			 tensor_declare(4,yyvsp[-6].string,yyvsp[-1].string);;
    break;}
case 22:
#line 155 "grpp_yacc.y"
{tensor_declare(1,yyvsp[-8].string,yyvsp[-7].string);
			 tensor_declare(1,yyvsp[-8].string,yyvsp[-5].string);
			 tensor_declare(1,yyvsp[-8].string,yyvsp[-3].string);
			 tensor_declare(1,yyvsp[-8].string,yyvsp[-1].string);;
    break;}
case 23:
#line 160 "grpp_yacc.y"
{tensor_declare(2,yyvsp[-8].string,yyvsp[-7].string);
			 tensor_declare(2,yyvsp[-8].string,yyvsp[-5].string);
			 tensor_declare(2,yyvsp[-8].string,yyvsp[-3].string);
			 tensor_declare(2,yyvsp[-8].string,yyvsp[-1].string);;
    break;}
case 24:
#line 165 "grpp_yacc.y"
{tensor_declare(3,yyvsp[-8].string,yyvsp[-7].string);
			 tensor_declare(3,yyvsp[-8].string,yyvsp[-5].string);
			 tensor_declare(3,yyvsp[-8].string,yyvsp[-3].string);
			 tensor_declare(3,yyvsp[-8].string,yyvsp[-1].string);;
    break;}
case 25:
#line 170 "grpp_yacc.y"
{tensor_declare(4,yyvsp[-8].string,yyvsp[-7].string);
			 tensor_declare(4,yyvsp[-8].string,yyvsp[-5].string);
			 tensor_declare(4,yyvsp[-8].string,yyvsp[-3].string);
			 tensor_declare(4,yyvsp[-8].string,yyvsp[-1].string);;
    break;}
case 26:
#line 175 "grpp_yacc.y"
{tensor_declare(1,yyvsp[-10].string,yyvsp[-9].string);
			 tensor_declare(1,yyvsp[-10].string,yyvsp[-7].string);
			 tensor_declare(1,yyvsp[-10].string,yyvsp[-5].string);
			 tensor_declare(1,yyvsp[-10].string,yyvsp[-3].string);
			 tensor_declare(1,yyvsp[-10].string,yyvsp[-1].string);;
    break;}
case 27:
#line 181 "grpp_yacc.y"
{tensor_declare(2,yyvsp[-10].string,yyvsp[-9].string);
			 tensor_declare(2,yyvsp[-10].string,yyvsp[-7].string);
			 tensor_declare(2,yyvsp[-10].string,yyvsp[-5].string);
			 tensor_declare(2,yyvsp[-10].string,yyvsp[-3].string);
			 tensor_declare(2,yyvsp[-10].string,yyvsp[-1].string);;
    break;}
case 28:
#line 187 "grpp_yacc.y"
{tensor_declare(3,yyvsp[-10].string,yyvsp[-9].string);
			 tensor_declare(3,yyvsp[-10].string,yyvsp[-7].string);
			 tensor_declare(3,yyvsp[-10].string,yyvsp[-5].string);
			 tensor_declare(3,yyvsp[-10].string,yyvsp[-3].string);
			 tensor_declare(3,yyvsp[-10].string,yyvsp[-1].string);;
    break;}
case 29:
#line 193 "grpp_yacc.y"
{tensor_declare(4,yyvsp[-10].string,yyvsp[-9].string);
			 tensor_declare(4,yyvsp[-10].string,yyvsp[-7].string);
			 tensor_declare(4,yyvsp[-10].string,yyvsp[-5].string);
			 tensor_declare(4,yyvsp[-10].string,yyvsp[-3].string);
			 tensor_declare(4,yyvsp[-10].string,yyvsp[-1].string);;
    break;}
case 30:
#line 199 "grpp_yacc.y"
{tensor_declare(1,yyvsp[-12].string,yyvsp[-11].string);
			 tensor_declare(1,yyvsp[-12].string,yyvsp[-9].string);
			 tensor_declare(1,yyvsp[-12].string,yyvsp[-7].string);
			 tensor_declare(1,yyvsp[-12].string,yyvsp[-5].string);
			 tensor_declare(1,yyvsp[-12].string,yyvsp[-3].string);
			 tensor_declare(1,yyvsp[-12].string,yyvsp[-1].string);;
    break;}
case 31:
#line 206 "grpp_yacc.y"
{tensor_declare(2,yyvsp[-12].string,yyvsp[-11].string);
			 tensor_declare(2,yyvsp[-12].string,yyvsp[-9].string);
			 tensor_declare(2,yyvsp[-12].string,yyvsp[-7].string);
			 tensor_declare(2,yyvsp[-12].string,yyvsp[-5].string);
			 tensor_declare(2,yyvsp[-12].string,yyvsp[-3].string);
			 tensor_declare(2,yyvsp[-12].string,yyvsp[-1].string);;
    break;}
case 32:
#line 213 "grpp_yacc.y"
{tensor_declare(3,yyvsp[-12].string,yyvsp[-11].string);
			 tensor_declare(3,yyvsp[-12].string,yyvsp[-9].string);
			 tensor_declare(3,yyvsp[-12].string,yyvsp[-7].string);
			 tensor_declare(3,yyvsp[-12].string,yyvsp[-5].string);
			 tensor_declare(3,yyvsp[-12].string,yyvsp[-3].string);
			 tensor_declare(3,yyvsp[-12].string,yyvsp[-1].string);;
    break;}
case 33:
#line 220 "grpp_yacc.y"
{tensor_declare(4,yyvsp[-12].string,yyvsp[-11].string);
			 tensor_declare(4,yyvsp[-12].string,yyvsp[-9].string);
			 tensor_declare(4,yyvsp[-12].string,yyvsp[-7].string);
			 tensor_declare(4,yyvsp[-12].string,yyvsp[-5].string);
			 tensor_declare(4,yyvsp[-12].string,yyvsp[-3].string);
			 tensor_declare(4,yyvsp[-12].string,yyvsp[-1].string);;
    break;}
case 34:
#line 229 "grpp_yacc.y"
{total_antisymmetric_tensor(yyvsp[-2].string);;
    break;}
case 35:
#line 231 "grpp_yacc.y"
{printf(" %s;\n",yyvsp[-1].string); no_lines_out++;;
    break;}
case 36:
#line 233 "grpp_yacc.y"
{vector_equation(yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 37:
#line 235 "grpp_yacc.y"
{equation_vector(yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 38:
#line 237 "grpp_yacc.y"
{tensor_eq_array(1,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 39:
#line 239 "grpp_yacc.y"
{tensor_eq_array(2,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 40:
#line 241 "grpp_yacc.y"
{tensor_eq_array(3,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 41:
#line 243 "grpp_yacc.y"
{tensor_eq_array(4,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 42:
#line 245 "grpp_yacc.y"
{array_eq_expression(yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 43:
#line 247 "grpp_yacc.y"
{tensor_eq_expression(0,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 44:
#line 249 "grpp_yacc.y"
{tensor_eq_expression(1,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 45:
#line 251 "grpp_yacc.y"
{tensor_eq_expression(2,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 46:
#line 253 "grpp_yacc.y"
{tensor_eq_expression(3,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 47:
#line 255 "grpp_yacc.y"
{tensor_eq_expression(4,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 48:
#line 257 "grpp_yacc.y"
{symten_eq_expression(0,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 49:
#line 259 "grpp_yacc.y"
{symten_eq_expression(1,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 50:
#line 261 "grpp_yacc.y"
{symten_eq_expression(2,yyvsp[-3].string,yyvsp[-2].string,yyvsp[-1].string);;
    break;}
case 51:
#line 265 "grpp_yacc.y"
{
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 52:
#line 269 "grpp_yacc.y"
{
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 53:
#line 273 "grpp_yacc.y"
{
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 54:
#line 277 "grpp_yacc.y"
{
			 expand_dummy_indices(1,yyvsp[0].string);
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 55:
#line 282 "grpp_yacc.y"
{
			 expand_dummy_indices(2,yyvsp[0].string);
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 56:
#line 287 "grpp_yacc.y"
{
			 expand_dummy_indices(3,yyvsp[0].string);
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 57:
#line 292 "grpp_yacc.y"
{
			 expand_dummy_indices(4,yyvsp[0].string);
			 strcpy(yyval.string,yyvsp[0].string);
			;
    break;}
case 58:
#line 297 "grpp_yacc.y"
{ covariant_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 59:
#line 299 "grpp_yacc.y"
{ covariant_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 60:
#line 301 "grpp_yacc.y"
{ covariant_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 61:
#line 303 "grpp_yacc.y"
{ covariant_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 62:
#line 305 "grpp_yacc.y"
{ lie_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 63:
#line 307 "grpp_yacc.y"
{ lie_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 64:
#line 309 "grpp_yacc.y"
{ lie_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 65:
#line 311 "grpp_yacc.y"
{ lie_derivative(yyvsp[-3].string,yyvsp[-1].string,yyval.string); ;
    break;}
case 66:
#line 313 "grpp_yacc.y"
{
			 if (!sum_expressions(yyval.string, yyvsp[-2].string, '+', yyvsp[0].string))
			   {
			    strcpy(msg,"Tensor Summation Error");
			    strcat(msg,yyvsp[-2].string);
			    strcat(msg," + ");
			    strcat(msg,yyvsp[0].string);
			    yyerror(msg);
			   }
			;
    break;}
case 67:
#line 324 "grpp_yacc.y"
{
			 if (!sum_expressions(yyval.string, yyvsp[-2].string, '-', yyvsp[0].string))
			   {
			    strcpy(msg,"Tensor Subtraction Error");
			    strcat(msg,yyvsp[-2].string);
			    strcat(msg," - ");
			    strcat(msg,yyvsp[0].string);
			    yyerror(msg);
			   }
			;
    break;}
case 68:
#line 335 "grpp_yacc.y"
{
			 if (!tensor_product(yyval.string,yyvsp[-2].string,yyvsp[0].string))
			   {
			    strcpy(msg,"Tensor Product Error");
			    strcat(msg,yyvsp[-2].string);
			    strcat(msg," * ");
			    strcat(msg,yyvsp[0].string);
			    yyerror(msg);
			   }
			;
    break;}
case 69:
#line 346 "grpp_yacc.y"
{
			    strcpy(yyval.string,yyvsp[-2].string);
			    strcat(yyval.string," / ");
			    strcat(yyval.string,yyvsp[0].string);
			;
    break;}
case 70:
#line 352 "grpp_yacc.y"
{
			 strcpy(yyval.string,"-");
			 strcat(yyval.string,yyvsp[0].string);
			;
    break;}
case 71:
#line 357 "grpp_yacc.y"
{
			 if ((yyvsp[-1].string[0] == '(')&&(yyvsp[-1].string[strlen(yyvsp[-1].string)-1] == ')'))
			   strcpy(yyval.string,yyvsp[-1].string);
			 else
			   {
			    strcpy(yyval.string,"( ");
			    strcat(yyval.string,yyvsp[-1].string);
			    strcat(yyval.string," )");
			   }
			;
    break;}
}
   /* the action file gets copied in in place of this dollarsign */
#line 543 "/usr/lib/bison.simple"

  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YYLSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = yylloc.first_line;
      yylsp->first_column = yylloc.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;

yyerrlab:   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      yyerror(msg);
	      free(msg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror("parse error");
    }

  goto yyerrlab1;
yyerrlab1:   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Discarding token %d (%s).\n", yychar, yytname[yychar1]);
#endif

      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;

yyerrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) goto yydefault;
#endif

yyerrpop:   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

yyerrhandle:

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;

 yyacceptlab:
  /* YYACCEPT comes here.  */
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
  return 0;

 yyabortlab:
  /* YYABORT comes here.  */
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
  return 1;
}
#line 369 "grpp_yacc.y"


/* Build the typedef tensor structures */

void tensor_structures(void)
{
 FILE *tf;
 int i,j,k,l,ii,jj,kk,ll;
 int rank,coord;
 char chr_i,chr_j,chr_k,chr_l;
 char element[5] = "";
 char *which_i_coord,*which_j_coord,*which_k_coord,*which_l_coord;
 static char prespace[] = "               ";
 
 if (coord_flag && index_flag) {
 if ((tf = fopen(hdrfile,"w")) == NULL)
   {
    fprintf(stderr,"Error: unable to create grpptd.h\n");
    exit(1);
   }
 else {

 /* Write out a header comment */
 fprintf(tf,"/* General Relativity PreProcessor */\n");
 fprintf(tf,"/*         Version 2.1             */\n");
 fprintf(tf,"/*                                 */\n");
 fprintf(tf,"/*     Copyright (c) 1992-2000     */\n");
 fprintf(tf,"/*     by James Kent Blackburn     */\n");
 fprintf(tf,"/*     All Rights Reserved         */\n");
 fprintf(tf,"/*                                 */\n");
 fprintf(tf,"/* GRPP header file: %s */\n\n",hdrfile);
 fprintf(tf,"\n#define GRPP_DIMENSION\t\t %d \n",n_dim);
 fprintf(tf,"#define GRPP_INDICES\t\t \"%s\" \n",index_list);
 fprintf(tf,"#define GRPP_COORDINATES\t \"%s\" \n\n",coord_list);
 fprintf(tf,"\n#define COMPONENT\t\t double \n");

 /* Write out the Rank 1 Tensor Structures */
 element[1] = '\0';
 for (ii=0; ii < 2; ii++)
    {
     if (ii == 0)
       {which_i_coord = lower_coord; chr_i = 'c';}
     else
       {which_i_coord = upper_coord; chr_i = 'C';}
     fprintf(tf,"\ntypedef struct { \n");
     fprintf(tf,"%s COMPONENT ",prespace);
     for (i=0; i<n_dim; i++)
        {
         element[0] = which_i_coord[i];
         if (i < (n_dim-1))
           fprintf(tf,"%s,",element);
         else
           fprintf(tf,"%s;\n",element);
        }
     fprintf(tf,"%s} TENSOR_%c;\n",prespace,chr_i);
    }

 /* Write out the Rank 2 Tensor Structures */
 element[2] = '\0';
 for (ii=0; ii < 2; ii++) {
 for (jj=0; jj < 2; jj++)
    {
     if (ii == 0)
       {which_i_coord = lower_coord; chr_i = 'c';}
     else
       {which_i_coord = upper_coord; chr_i = 'C';}
     if (jj == 0)
       {which_j_coord = lower_coord; chr_j = 'c';}
     else
       {which_j_coord = upper_coord; chr_j = 'C';}
     fprintf(tf,"\ntypedef struct { \n");
     for (i=0; i<n_dim; i++) {
     fprintf(tf,"%s COMPONENT ",prespace);
     for (j=0; j<n_dim; j++)
        {
         element[0] = which_i_coord[i];
         element[1] = which_j_coord[j];
         if (j < (n_dim-1))
           fprintf(tf,"%s,",element);
         else
           fprintf(tf,"%s;\n",element);
        } }
     fprintf(tf,"%s} TENSOR_%c%c;\n",prespace,chr_i,chr_j);
    } }

 /* Write out the Rank 3 Tensor Structures */
 element[3] = '\0';
 for (ii=0; ii < 2; ii++) {
 for (jj=0; jj < 2; jj++) {
 for (kk=0; kk < 2; kk++)
    {
     if (ii == 0)
       {which_i_coord = lower_coord; chr_i = 'c';}
     else
       {which_i_coord = upper_coord; chr_i = 'C';}
     if (jj == 0)
       {which_j_coord = lower_coord; chr_j = 'c';}
     else
       {which_j_coord = upper_coord; chr_j = 'C';}
     if (kk == 0)
       {which_k_coord = lower_coord; chr_k = 'c';}
     else
       {which_k_coord = upper_coord; chr_k = 'C';}
     fprintf(tf,"\ntypedef struct { \n");
     for (i=0; i<n_dim; i++) {
     for (j=0; j<n_dim; j++) {
     fprintf(tf,"%s COMPONENT ",prespace);
     for (k=0; k<n_dim; k++)
        {
         element[0] = which_i_coord[i];
         element[1] = which_j_coord[j];
         element[2] = which_k_coord[k];
         if (k < (n_dim-1))
           fprintf(tf,"%s,",element);
         else
           fprintf(tf,"%s;\n",element);
        } } }
     fprintf(tf,"%s} TENSOR_%c%c%c;\n",prespace,chr_i,chr_j,chr_k);
    } } }

 /* Write out the Rank 4 Tensor Structures */
 element[4] = '\0';
 for (ii=0; ii < 2; ii++) {
 for (jj=0; jj < 2; jj++) {
 for (kk=0; kk < 2; kk++) {
 for (ll=0; ll < 2; ll++)
    {
     if (ii == 0)
       {which_i_coord = lower_coord; chr_i = 'c';}
     else
       {which_i_coord = upper_coord; chr_i = 'C';}
     if (jj == 0)
       {which_j_coord = lower_coord; chr_j = 'c';}
     else
       {which_j_coord = upper_coord; chr_j = 'C';}
     if (kk == 0)
       {which_k_coord = lower_coord; chr_k = 'c';}
     else
       {which_k_coord = upper_coord; chr_k = 'C';}
     if (ll == 0)
       {which_l_coord = lower_coord; chr_l = 'c';}
     else
       {which_l_coord = upper_coord; chr_l = 'C';}
     fprintf(tf,"\ntypedef struct { \n");
     for (i=0; i<n_dim; i++) {
     for (j=0; j<n_dim; j++) {
     for (k=0; k<n_dim; k++) {
     fprintf(tf,"%s COMPONENT ",prespace);
     for (l=0; l<n_dim; l++)
        {
         element[0] = which_i_coord[i];
         element[1] = which_j_coord[j];
         element[2] = which_k_coord[k];
         element[3] = which_l_coord[l];
         if (l < (n_dim-1))
           fprintf(tf,"%s,",element);
         else
           fprintf(tf,"%s;\n",element);
        } } } }
     fprintf(tf,"%s} TENSOR_%c%c%c%c;\n",prespace,chr_i,chr_j,chr_k,chr_l);
    } } } }
 fprintf(tf,"\n#undef COMPONENT \n");
 fclose(tf);}}
}

/* Write out the tensor component declarations statements */

void tensor_declare(int rank, char *type, char *tensor)
{
 int	number_indices,i,j;
 char	chr_i,chr_j,chr_k,chr_l;
 char	t_name[32] = "";
 char   t_temp[32] = "";
 char	t_indx[6] = "";
 i = 0; j = 0;
 while (tensor[i] != '\0')
      {
       if ((tensor[i] != '(') && (tensor[i] != ')'))
	 {
	  t_temp[j] = tensor[i];
	  j++;
	 }
       i++;
      }
 t_temp[j] = '\0';
 number_indices = parse_tensor(rank,t_temp,t_name,t_indx);
 switch (rank)
       {
	case 1:
	  if (isupper(t_indx[0])) chr_i = 'C';
	  else chr_i = 'c';
	  printf(" TENSOR_%c %s ;\n",chr_i,t_name);
	  break;
	case 2:
	  if (isupper(t_indx[0])) chr_i = 'C';
	  else chr_i = 'c';
	  if (isupper(t_indx[1])) chr_j = 'C';
	  else chr_j = 'c';
	  printf(" TENSOR_%c%c %s ;\n",chr_i,chr_j,t_name);
	  break;
	case 3:
	  if (isupper(t_indx[0])) chr_i = 'C';
	  else chr_i = 'c';
	  if (isupper(t_indx[1])) chr_j = 'C';
	  else chr_j = 'c';
	  if (isupper(t_indx[2])) chr_k = 'C';
	  else chr_k = 'c';
 	  printf(" TENSOR_%c%c%c %s ;\n",chr_i,chr_j,chr_k,t_name);
	  break;
	case 4:
	  if (isupper(t_indx[0])) chr_i = 'C';
	  else chr_i = 'c';
	  if (isupper(t_indx[1])) chr_j = 'C';
	  else chr_j = 'c';
	  if (isupper(t_indx[2])) chr_k = 'C';
	  else chr_k = 'c';
	  if (isupper(t_indx[3])) chr_l = 'C';
	  else chr_l = 'c';
	  printf(" TENSOR_%c%c%c%c %s ;\n",chr_i,chr_j,chr_k,chr_l,t_name);
	  break;
       }
 no_lines_out++;
}

/* Build the totally antisymmetric tensor */

void total_antisymmetric_tensor(char *prefix)
{
 int	i,j,k,l,len,ndollar,value,list[5];
 char	epsilon[32] = "";
 char   msg[80] = "";
 char   svalue[3] = "";

 strcpy(epsilon,prefix);
 switch (n_dim)
       {
	case 1:
	   strcat(epsilon,"_C");
	break;
	case 2:
	   strcat(epsilon,"_CC");
	break;
	case 3:
	   strcat(epsilon,"_CCC");
	break;
	case 4:
	   strcat(epsilon,"_CCCC");
	break;
       }
 len = strlen(epsilon);
 epsilon[len] = '$';
 ndollar = len;
 len++;
 epsilon[len] = '\0';
 len++;
 switch (n_dim)
       {
	case 1:
	  strcpy(msg,"Rank 1 Totally Antisymmetric Tensor Error");
	  yyerror(msg);
	break ;
	case 2:
	  epsilon[ndollar+1] = upper_index[0];
	  epsilon[ndollar+2] = upper_index[1];
	  epsilon[ndollar+3] = '\0';
          tensor_eq_expression(2,epsilon,"="," 0");
	  epsilon[ndollar] = '.';
	  for (i=0;i<n_dim;i++)
	     for (j=0;j<n_dim;j++)
		{
		 list[0] = i; list[1] = j;
		 if ( !(i==j) )
		   {
		    value = permutation_sort(n_dim,list);
		    if ( value == -1 )
		      strcpy(svalue,"-1");
		    else
		      strcpy(svalue," 1");
		    epsilon[ndollar+1] = upper_coord[i];
		    epsilon[ndollar+2] = upper_coord[j];
		    epsilon[ndollar+3] = '\0';
		    tensor_eq_expression(0,epsilon,"=",svalue);
		   }
		}
	break;
	case 3:
	  epsilon[ndollar+1] = upper_index[0];
	  epsilon[ndollar+2] = upper_index[1];
	  epsilon[ndollar+3] = upper_index[2];
	  epsilon[ndollar+4] = '\0';
          tensor_eq_expression(3,epsilon,"="," 0");
	  epsilon[ndollar] = '.';
	  for (i=0;i<n_dim;i++)
	     for (j=0;j<n_dim;j++)
		for (k=0;k<n_dim;k++)
		   {
		    list[0] = i; list[1] = j; list[2] = k;
		    if ( !((i==j)||(i==k)||(j==k)) )
		      {
		       value = permutation_sort(n_dim,list);
		       if ( value == -1 )
			 strcpy(svalue,"-1");
		       else
			 strcpy(svalue," 1");
		       epsilon[ndollar+1] = upper_coord[i];
		       epsilon[ndollar+2] = upper_coord[j];
		       epsilon[ndollar+3] = upper_coord[k];
		       epsilon[ndollar+4] = '\0';
		       tensor_eq_expression(0,epsilon,"=",svalue);
		      }
		   }
	break;
	case 4:
	  strcat(epsilon,upper_coord);
          tensor_eq_expression(4,epsilon,"="," 0");
	  epsilon[ndollar] = '.';
	  for (i=0;i<n_dim;i++)
	     for (j=0;j<n_dim;j++)
		for (k=0;k<n_dim;k++)
		   for (l=0;l<n_dim;l++)
		      {
		       list[0] = i; list[1] = j; list[2] = k; list[3] = l;
		       if ( !((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) )
		         {
		          value = permutation_sort(n_dim,list);
		          if ( value == -1 )
			    strcpy(svalue,"-1");
		          else
			    strcpy(svalue," 1");
		          epsilon[ndollar+1] = upper_coord[i];
		          epsilon[ndollar+2] = upper_coord[j];
		          epsilon[ndollar+3] = upper_coord[k];
		          epsilon[ndollar+4] = upper_coord[l];
		          epsilon[ndollar+5] = '\0';
		          tensor_eq_expression(0,epsilon,"=",svalue);
		      }
		   }
	break;
       }
}

/* Use bubble sort to determine if this is a even or odd permutation */

int permutation_sort(int elements, int *list)
{
 int	i,j,holder,permutations;

 permutations = 0;
 for (i=0;i<elements-1;i++)
    {
     for (j=i+1;j<elements;j++)
	{
	 if ( list[i] > list[j] )
	   {
	    holder = list[i];
	    list[i] = list[j];
	    list[j] = holder;
	    permutations++;
	   }
	}
    }
 if (permutations % 2)
   return(-1);
 else
  return(1);
}

/* Make a vector assignment to a rank 1 tensor*/

void vector_equation(char *tensor, char *eq, char *vector)
{
 int	i,j,m,len,mrank;
 char	*dollar,actual_ind,subst;
 char   name[32] = "";
 char   vec[160] = "";
 char	*component,
	delimit[]=",",
	msg[] = "Vector Component Error - check dimension";

 if (gr_debug) fprintf(stderr,"TRACE: %s %s %s\n",tensor,eq,vector);
 for (j=0;j<32;j++)
    name[j] = '\0';
 for (j=0;j<160;j++)
    vec[j] = '\0';
 len = strlen(vector) - 2;
 strncpy(vec,&vector[1],len);
 mrank = is_index(tensor) + is_coord(tensor);
 strcpy(name,tensor);
 dollar = (char *)index(name,'$');
 *dollar = '.';
 len = m = 0;
 do {
     if (index(index_list,*(dollar+m+1))) len = m+1;
     m++;
    } while (len==0);

 component = strtok(vec,delimit);
 actual_ind = *(dollar+len);
 for (i=0;i<n_dim;i++)
    {
     if (islower(actual_ind))
       subst = lower_coord[i];
     else
       subst = upper_coord[i];
     *(dollar+len) = subst;
     if ( component != NULL)
       print_equation(name,eq,component);
     else
       yyerror(msg);
     component = strtok(NULL,delimit);
    }
}

/* Make a vector assignment to a rank 1 tensor*/

void equation_vector(char *vector, char *eq, char *tensor)
{
 int	i,j,m,len,mrank;
 char	*dollar,actual_ind,subst;
 char   name[32] = "";
 char   vec[160] = "";
 char	*component,
	delimit[]=",",
	msg[] = "Vector Component Error - check dimension";

 if (gr_debug) fprintf(stderr,"TRACE: %s %s %s\n",tensor,eq,vector);
 for (j=0;j<32;j++)
    name[j] = '\0';
 for (j=0;j<160;j++)
    vec[j] = '\0';
 len = strlen(vector) - 2;
 strncpy(vec,&vector[1],len);
 mrank = is_index(tensor) + is_coord(tensor);
 strcpy(name,tensor);
 dollar = (char *)index(name,'$');
 *dollar = '.';
 len = m = 0;
 do {
     if (index(index_list,*(dollar+m+1))) len = m+1;
     m++;
    } while (len==0);

 component = strtok(vec,delimit);
 actual_ind = *(dollar+len);
 for (i=0;i<n_dim;i++)
    {
     if (islower(actual_ind))
       subst = lower_coord[i];
     else
       subst = upper_coord[i];
     *(dollar+len) = subst;
     if ( component != NULL)
       print_equation(component,eq,name);
     else
       yyerror(msg);
     component = strtok(NULL,delimit);
    }
}

/* Assign array components to a tensor */

void tensor_eq_array(int rank, char *tensor, char *eq, char *array)
{
 int	i,j,k,l,m,n,len,number_indices,mrank;
 char	*dollar,actual_ind,subst;
 char   temp[32] = "";
 char   tend[32] = "";
 char	tensor_name[32] = "";
 char   array_name[32] = "";
 char   t_indx[25] = "";
 char	tindex[5] = "";
 char   buffer[8192] = "";

 if (gr_debug) fprintf(stderr,"TRACE: %s %s %s\n",tensor,eq,array);
 for (j=0;j<25;j++) 
    t_indx[j] = '\0';
 for (j=0;j<32;j++)
    {
     temp[j] = '\0'; 
     tensor_name[j] = '\0';
     array_name[j] = '\0';
    }
 len = 0;
 if (rank > 0)
   {
    mrank = is_index(tensor) + is_coord(tensor);
    dollar = (char *)index(tensor,'$');
    strcpy(tend,(dollar+1));
    n = 0;
    for (m=0; m<mrank; m++)
       if (index(index_list, *(dollar+m+1)))
	 {tindex[n] = *(dollar+m+1); n++;}
    number_indices = parse_tensor(rank,tensor,tensor_name,t_indx);
    strcpy(temp,tensor_name);
    len = strlen(temp);
    temp[len] = '.';
    len++;
    temp[len] = '\0';
    strcat(temp,tend);
   }
 else
    number_indices = 0;
 switch (number_indices)
       {
	case 1:
	  for (i=0;i<n_dim;i++)
	     {
	      for (m=0;m<mrank;m++)
		 {
		  actual_ind = *(dollar+m+1);
		  if (index(index_list,actual_ind))
		    {
		     if (islower(actual_ind))
		       subst = lower_coord[i];
		     else
		       subst = upper_coord[i];
		     temp[len+m] = subst;
		    }
		 }
	      array_elements(array_name,array,1,i,0,0,0);
	      print_equation(temp,eq,array_name);
	     }
	  break;
	case 2:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (m=0;m<mrank;m++)
		     {
		      actual_ind = *(dollar+m+1);
		      if (index(index_list,actual_ind))
			{
		         if (isupper(actual_ind))
			   {
			    if (actual_ind == toupper(tindex[0]))
			      subst = upper_coord[i];
			    else
			      subst = upper_coord[j];
			   }
		         else
			   {
			    if (actual_ind == tolower(tindex[0]))
			      subst = lower_coord[i];
			    else
			      subst = lower_coord[j];
			   }
		         temp[len+m] = subst;
			}
		     }
	          array_elements(array_name,array,2,i,j,0,0);
		  print_equation(temp,eq,array_name);
		 }
	     }
 	  break;
	case 3:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (k=0;k<n_dim;k++)
		     {
		      for (m=0;m<mrank;m++)
			 {
			  actual_ind = *(dollar+m+1);
			  if (index(index_list,actual_ind))
			    {
			     if (isupper(actual_ind))
			       {
			        if (actual_ind == toupper(tindex[0]))
			          subst = upper_coord[i];
			        else if (actual_ind == toupper(tindex[1]))
			          subst = upper_coord[j];
			        else
			          subst = upper_coord[k];
			       }
		             else
			       {
			        if (actual_ind == tolower(tindex[0]))
			          subst = lower_coord[i];
			        else if (actual_ind == tolower(tindex[1]))
			          subst = lower_coord[j];
			        else
			          subst = lower_coord[k];
			       }
		             temp[len+m] = subst;
			    }
		         }
	              array_elements(array_name,array,3,i,j,k,0);
		      print_equation(temp,eq,array_name);
		     }
		 }
	     }
	  break;
	case 4:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (k=0;k<n_dim;k++)
		     {
		      for (l=0;l<n_dim;l++)
			 {
		          for (m=0;m<mrank;m++)
			     {
			      actual_ind = *(dollar+m+1);
			      if (index(index_list,actual_ind)){
			      if (isupper(actual_ind))
			        {
			         if (actual_ind == toupper(tindex[0]))
			           subst = upper_coord[i];
			         else if (actual_ind == toupper(tindex[1]))
			           subst = upper_coord[j];
			         else if (actual_ind == toupper(tindex[2]))
			           subst = upper_coord[k];
			         else
			           subst = upper_coord[l];
			        }
		              else
			        {
			         if (actual_ind == tolower(tindex[0]))
			           subst = lower_coord[i];
			         else if (actual_ind == tolower(tindex[1]))
			           subst = lower_coord[j];
			         else if (actual_ind == tolower(tindex[2]))
			           subst = lower_coord[k];
			         else
			           subst = lower_coord[l];
			        }
		              temp[len+m] = subst;}
		             }
	                  array_elements(array_name,array,4,i,j,k,l);
	      		  print_equation(temp,eq,array_name);
			 }
		     }
		 }
	     }
	  break;
       }
}

/* Assign expression to an array */

void array_eq_expression(char *array, char *eq, char *express)
{
 int  i,j,k,l,m,n,len,number_indices;
 char *dollar,actual_ind,subst;
 char	array_name[32] = "";
 char   free[12] = "";
 char	buffer[8192] = "";

 if (gr_debug) fprintf(stderr,"TRACE: %s %s %s\n",array,eq,express);
 for (j=0;j<32;j++)
    array_name[j] = '\0';
 len = 0;
 number_indices = get_free_indices(express, free);
 dollar = &free[0];
 switch (number_indices)
       {
	case 1:
	  for (i=0;i<n_dim;i++)
	     {
	      strncpy(buffer,express,8192);
	      for (m=0;m<number_indices;m++)
		 {
		  actual_ind = *(dollar+m);
		  if (islower(actual_ind))
		    subst = lower_coord[i];
		  else
		    subst = upper_coord[i];
	          expression_replace(actual_ind,subst,buffer);
		 }
	      array_elements(array_name,array,1,i,0,0,0);
	      print_equation(array_name,eq,buffer);
	     }
	  break;
	case 2:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  strncpy(buffer,express,8192);
		  for (m=0;m<number_indices;m++)
		     {
		      actual_ind = *(dollar+m);
		      if (isupper(actual_ind))
			{
			 if (actual_ind == toupper(free[0]))
			   subst = upper_coord[i];
			 else
			   subst = upper_coord[j];
			}
		      else
			{
			 if (actual_ind == tolower(free[0]))
			   subst = lower_coord[i];
			 else
			   subst = lower_coord[j];
			}
	              expression_replace(actual_ind,subst,buffer);
		     }
		  array_elements(array_name,array,2,i,j,0,0);
		  print_equation(array_name,eq,buffer);
		 }
	     }
 	  break;
	case 3:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (k=0;k<n_dim;k++)
		     {
		      strncpy(buffer,express,8192);
		      for (m=0;m<number_indices;m++)
			 {
			  actual_ind = *(dollar+m);
			  if (isupper(actual_ind))
			    {
			     if (actual_ind == toupper(free[0]))
			       subst = upper_coord[i];
			     else if (actual_ind == toupper(free[1]))
			       subst = upper_coord[j];
			     else
			       subst = upper_coord[k];
			    }
		          else
			    {
			     if (actual_ind == tolower(free[0]))
			       subst = lower_coord[i];
			     else if (actual_ind == tolower(free[1]))
			       subst = lower_coord[j];
			     else
			       subst = lower_coord[k];
			    }
	                  expression_replace(actual_ind,subst,buffer);
		         }
		      array_elements(array_name,array,3,i,j,k,0);
		      print_equation(array_name,eq,buffer);
		     }
		 }
	     }
	  break;
	case 4:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (k=0;k<n_dim;k++)
		     {
		      for (l=0;l<n_dim;l++)
			 {
			  strncpy(buffer,express,8192);
		          for (m=0;m<number_indices;m++)
			     {
			      actual_ind = *(dollar+m);
			      if (isupper(actual_ind))
			        {
			         if (actual_ind == toupper(free[0]))
			           subst = upper_coord[i];
			         else if (actual_ind == toupper(free[1]))
			           subst = upper_coord[j];
			         else if (actual_ind == toupper(free[2]))
			           subst = upper_coord[k];
			         else
			           subst = upper_coord[l];
			        }
		              else
			        {
			         if (actual_ind == tolower(free[0]))
			           subst = lower_coord[i];
			         else if (actual_ind == tolower(free[1]))
			           subst = lower_coord[j];
			         else if (actual_ind == tolower(free[2]))
			           subst = lower_coord[k];
			         else
			           subst = lower_coord[l];
			        }
	                      expression_replace(actual_ind,subst,buffer);
		             }
		          array_elements(array_name,array,4,i,j,k,l);
	      		  print_equation(array_name,eq,buffer);
			 }
		     }
		 }
	     }
	  break;
       }
}

/* build the array name for a particular element */

void array_elements(char *elem, char *arr, int rk, int i, int j, int k, int l)
{
 int m = 0;

 do {
     elem[m] = arr[m];
     m++;
    }
 while ( arr[m] != '[');
 switch (rk)
       {
	case 1:
	  elem[m] = '[';
	  elem[m+1] = (char)('0' + i);
	  elem[m+2] = ']';
	  m += 3;
	  break;
	case 2:
	  elem[m] = '[';
	  elem[m+1] = (char)('0' + i);
	  elem[m+2] = ']';
	  elem[m+3] = '[';
	  elem[m+4] = (char)('0' + j);
	  elem[m+5] = ']';
	  m += 6;
	  break;
	case 3:
	  elem[m] = '[';
	  elem[m+1] = (char)('0' + i);
	  elem[m+2] = ']';
	  elem[m+3] = '[';
	  elem[m+4] = (char)('0' + j);
	  elem[m+5] = ']';
	  elem[m+6] = '[';
	  elem[m+7] = (char)('0' + k);
	  elem[m+8] = ']';
	  m += 9;
	  break;
	case 4:
	  elem[m] = '[';
	  elem[m+1] = (char)('0' + i);
	  elem[m+2] = ']';
	  elem[m+3] = '[';
	  elem[m+4] = (char)('0' + j);
	  elem[m+5] = ']';
	  elem[m+6] = '[';
	  elem[m+7] = (char)('0' + k);
	  elem[m+8] = ']';
	  elem[m+9] = '[';
	  elem[m+10] = (char)('0' + l);
	  elem[m+11] = ']';
	  m += 12;
	  break;
       }
 elem[m] = '\0';
}

/* Sum two expressions */

int sum_expressions(char *exp, char *term1, char sign, char *term2)
{
 int rank1,rank2;
 
 if (gr_debug) fprintf(stderr,"TRACE: %s %c %s ->\n",term1,sign,term2);
 rank1 = get_rank(term1);
 rank2 = get_rank(term2);
 if ( rank1 == rank2)
   {
    strcpy(exp,term1);
    if (sign == '+')
      {
       strcat(exp," + ");
       strcat(exp,term2);
      }
    else
      {
       strcat(exp," - ");
       strcat(exp,term2);
      }
    if (gr_debug) fprintf(stderr,"       %s\n",exp);
    return 1;
   }
 else
   return 0;
}

/* Expand the dummy indices in terms of the coordinates in a tensor */

void expand_dummy_indices(int rank, char *tensor)
{
 int	i,j,m,len,num_pairs,total;
 char	*dollar,actual_ind;
 char   dummy[25] = "";
 char   temp[32] = "";
 char   buffer[8192] = "";

 if (gr_debug) fprintf(stderr,"TRACE: %s ->\n",tensor);
 for (i=0;i<25;i++) dummy[i] = '\0';
 for (i=0;i<8192;i++) buffer[i] = '\0';
 num_pairs = get_dummy_indices(tensor,dummy) / 2;
 if (num_pairs > 0)
   {
    dollar = (char *)index(tensor,'$');
    len = strlen(tensor) - rank;
    total = 0;
    switch (num_pairs)
       {
	case 1:
	  for (i=0;i<n_dim;i++)
	     {
	      strcpy(temp,tensor);
	      for (m=0;m<rank;m++)
		 {
		  actual_ind = *(dollar+m+1);
		  if ((dummy[0] == actual_ind)||(dummy[1] == actual_ind))
		    {
		     if (islower(actual_ind))
		       temp[len+m] = lower_coord[i];
		     else
		       temp[len+m] = upper_coord[i];
		    }
		 }
	      strcat(buffer,temp);
	      if (i < (n_dim-1)) strcat(buffer," + ");
	     }
	 break;
	 case 2:
	   for (i=0;i<n_dim;i++)
	      {
	       for (j=0;j<n_dim;j++)
	          {
	           strcpy(temp,tensor);
	           for (m=0;m<rank;m++)
		      {
		       actual_ind = *(dollar+m+1);
		       if ((dummy[0]==actual_ind)||(dummy[1]==actual_ind))
		         {
		          if (islower(actual_ind))
		            temp[len+m] = lower_coord[i];
		          else
		            temp[len+m] = upper_coord[i];
		         }
		       else if ((dummy[2]==actual_ind)||(dummy[3]==actual_ind))
			 {
		          if (islower(actual_ind))
		            temp[len+m] = lower_coord[j];
		          else
		            temp[len+m] = upper_coord[j];
		 	 }
		      }
		   total++;
	           strcat(buffer,temp);
	           if (total < (n_dim*n_dim)) strcat(buffer," + ");
		  }
	      }
	 break;
       }
    strcpy(tensor,"( ");
    strcat(tensor,buffer);
    strcat(tensor," )");
   }
 if (gr_debug) fprintf(stderr,"       %s\n",tensor);
}

/* Product two tensors, checking for both inner and outer products */

int tensor_product(char *exp, char *exp1, char *exp2)
{
 int	i,j,k,l,m,len1,len2,num_pairs,total,rank,rk1,rk2;
 char	*dol1,*dol2,subst1,subst2,act_ind;
 char   dummy[25] = "";
 char	free1[25] = "";
 char   free2[25] = "";
 char	psuedo[32] = "";
 char   buffer[8192] = "";
 char   temp1[8192] = "";
 char   temp2[8192] = "";

 for (i=0;i<8192;i++) 
    {
     buffer[i] = '\0';
     temp1[i] = '\0';
     temp2[i] = '\0';
     }
 for (i=0;i<25;i++)
    {
     dummy[i] = '\0';
     free1[i] = '\0';
     free2[i] = '\0';
    }
 for (i=0;i<32;i++) psuedo[i] = '\0';
 strcpy(psuedo,"psuedo$");
 strcpy(buffer,"( ");
 rk1 = get_free_indices(exp1, free1);
 rk2 = get_free_indices(exp2, free2);
 dol1 = &free1[0];
 dol2 = &free2[0];
 strcat(psuedo,free1);
 strcat(psuedo,free2);
 num_pairs = get_dummy_indices(psuedo,dummy) / 2;
 total = 0;
 switch (num_pairs)
    {
     case 0:
	strcat(buffer,exp1);
	strcat(buffer," * ");
	strcat(buffer,exp2);
     break;
     case 1:
	for (i=0;i<n_dim;i++)
	   {
	    strcpy(temp1,exp1);
	    strcpy(temp2,exp2);
	    for (m=0;m<rk1;m++)
	       {
		act_ind = *(dol1+m);
		if ((dummy[0] == act_ind)||(dummy[1] == act_ind))
		  {
		   if (islower(act_ind))
		     subst1 = lower_coord[i];
		   else
		     subst1 = upper_coord[i];
		  expression_replace(act_ind,subst1,temp1);
		  }
	       }
	    for (m=0;m<rk2;m++)
	       {
		act_ind = *(dol2+m);
		if ((dummy[0] == act_ind)||(dummy[1] == act_ind))
		  {
		   if (islower(act_ind))
		     subst2 = lower_coord[i];
		   else
		     subst2 = upper_coord[i];
		   expression_replace(act_ind,subst2,temp2);
		  }
	       }
	    strcat(buffer,temp1);
	    strcat(buffer," * ");
	    strcat(buffer,temp2);
	    if (i < (n_dim-1)) strcat(buffer," + ");
	   }
     break;
     case 2:
	for (i=0;i<n_dim;i++)
	   {
	    for (j=0;j<n_dim;j++)
	       {
		strcpy(temp1,exp1);
		strcpy(temp2,exp2);
	        for (m=0;m<rk1;m++)
		   {
		    act_ind = *(dol1+m);
		    if ((dummy[0]==act_ind)||(dummy[1]==act_ind))
		      {
		       if (islower(act_ind))
		         subst1 = lower_coord[i];
		       else
		         subst1 = upper_coord[i];
		       expression_replace(act_ind,subst1,temp1);
		      }
		    else if ((dummy[2]==act_ind)||(dummy[3]==act_ind))
		      {
		       if (islower(act_ind))
		         subst1 = lower_coord[j];
		       else
		         subst1 = upper_coord[j];
		       expression_replace(act_ind,subst1,temp1);
		      }
		   }
	        for (m=0;m<rk2;m++)
		   {
		    act_ind = *(dol2+m);
		    if ((dummy[0]==act_ind)||(dummy[1]==act_ind))
		      {
		       if (islower(act_ind))
		         subst2 = lower_coord[i];
		       else
		         subst2 = upper_coord[i];
		       expression_replace(act_ind,subst2,temp2);
		      }
		    else if ((dummy[2]==act_ind)||(dummy[3]==act_ind))
		      {
		       if (islower(act_ind))
		         subst2 = lower_coord[j];
		       else
		         subst2 = upper_coord[j];
		       expression_replace(act_ind,subst2,temp2);
		      }
		   }
	        strcat(buffer,temp1);
	        strcat(buffer," * ");
	        strcat(buffer,temp2);
	        total++;
	        if (total < (n_dim*n_dim)) strcat(buffer," + ");
	       }	
	   }
     break;
     case 3:
	for (i=0;i<n_dim;i++)
	   {
	    for (j=0;j<n_dim;j++)
	       {
		for (k=0;k<n_dim;k++)
		   {
		    strcpy(temp1,exp1);
		    strcpy(temp2,exp2);
	            for (m=0;m<rk1;m++)
		       {
		        act_ind = *(dol1+m);
		        if ((dummy[0]==act_ind)||(dummy[1]==act_ind))
		          {
		           if (islower(act_ind))
		             subst1 = lower_coord[i];
		           else
		             subst1 = upper_coord[i];
		           expression_replace(act_ind,subst1,temp1);
		          }
		        else if ((dummy[2]==act_ind)||(dummy[3]==act_ind))
		          {
		           if (islower(act_ind))
		             subst1 = lower_coord[j];
		           else
		             subst1 = upper_coord[j];
		           expression_replace(act_ind,subst1,temp1);
		 	  }
		        else 
		          {
		           if (islower(act_ind))
		             subst1 = lower_coord[k];
		           else
		             subst1 = upper_coord[k];
		           expression_replace(act_ind,subst1,temp1);
		 	  }
		       }
	            for (m=0;m<rk2;m++)
		       {
		        act_ind = *(dol2+m);
		        if ((dummy[0]==act_ind)||(dummy[1]==act_ind))
		          {
		           if (islower(act_ind))
		             subst2 = lower_coord[i];
		           else
		             subst2 = upper_coord[i];
		           expression_replace(act_ind,subst2,temp2);
		          }
		        else if ((dummy[2]==act_ind)||(dummy[3]==act_ind))
		          {
		           if (islower(act_ind))
		             subst2 = lower_coord[j];
		           else
		             subst2 = upper_coord[j];
		           expression_replace(act_ind,subst2,temp2);
		 	  }
		        else 
		          {
		           if (islower(act_ind))
		             subst2 = lower_coord[k];
		           else
		             subst2 = upper_coord[k];
		           expression_replace(act_ind,subst2,temp2);
		 	  }
		       }
	            strcat(buffer,temp1);
	            strcat(buffer," * ");
	            strcat(buffer,temp2);
	            total++;
	            if (total < (n_dim*n_dim*n_dim)) strcat(buffer," + ");
                   }
               }
           }
     break;
     case 4:
	for (i=0;i<n_dim;i++)
	   {
	    for (j=0;j<n_dim;j++)
	       {
		for (k=0;k<n_dim;k++)
		   {
		    for (l=0;l<n_dim;l++)
		       {
		        strcpy(temp1,exp1);
		        strcpy(temp2,exp2);
	                for (m=0;m<rk1;m++)
		           {
		            act_ind = *(dol1+m);
		            if ((dummy[0]==act_ind)||(dummy[1]==act_ind))
		              {
		               if (islower(act_ind))
		                 subst1 = lower_coord[i];
		               else
		                 subst1 = upper_coord[i];
		               expression_replace(act_ind,subst1,temp1);
		              }
		            else if ((dummy[2]==act_ind)||(dummy[3]==act_ind))
		              {
		               if (islower(act_ind))
		                 subst1 = lower_coord[j];
		               else
		                 subst1 = upper_coord[j];
		               expression_replace(act_ind,subst1,temp1);
		 	      }
		            else if ((dummy[4]==act_ind)||(dummy[5]==act_ind))
		              {
		               if (islower(act_ind))
		                 subst1 = lower_coord[k];
		               else
		                 subst1 = upper_coord[k];
		               expression_replace(act_ind,subst1,temp1);
		 	      }
		            else 
		              {
		               if (islower(act_ind))
		                 subst1 = lower_coord[l];
		               else
		                 subst1 = upper_coord[l];
		               expression_replace(act_ind,subst1,temp1);
		 	      }
		           }
	                for (m=0;m<rk2;m++)
		           {
		            act_ind = *(dol2+m);
		            if ((dummy[0]==act_ind)||(dummy[1]==act_ind))
		              {
		               if (islower(act_ind))
		                 subst2 = lower_coord[i];
		               else
		                 subst2 = upper_coord[i];
		               expression_replace(act_ind,subst2,temp2);
		              }
		            else if ((dummy[2]==act_ind)||(dummy[3]==act_ind))
		              {
		               if (islower(act_ind))
		                 subst2 = lower_coord[j];
		               else
		                 subst2 = upper_coord[j];
		               expression_replace(act_ind,subst2,temp2);
		 	      }
		            else if ((dummy[4]==act_ind)||(dummy[5]==act_ind))
		              {
		               if (islower(act_ind))
		                 subst2 = lower_coord[k];
		               else
		                 subst2 = upper_coord[k];
		               expression_replace(act_ind,subst2,temp2);
		 	      }
		            else 
		              {
		               if (islower(act_ind))
		                 subst2 = lower_coord[l];
		               else
		                 subst2 = upper_coord[l];
		               expression_replace(act_ind,subst2,temp2);
		 	      }
		           }
	                strcat(buffer,temp1);
	                strcat(buffer," * ");
	                strcat(buffer,temp2);
	                total++;
	                if (total<(n_dim*n_dim*n_dim*n_dim)) 
			  strcat(buffer," + ");
		       }
                   }
               }
           }
     break;
    }
 strcat(buffer," )");
 strcpy(exp,buffer);
 if (gr_debug)
   {
    fprintf(stderr,"TRACE: %s * %s ->\n",exp1,exp2);
    fprintf(stderr,"       %s\n",exp);
   }
 return 1;
}

/* Assign a symetric tensor to an expression component by component */

void symten_eq_expression(int rank, char *tensor, char *seq, char *expression)
{
 int  i,j,k,l,m,n,len,number_indices,mrank,loc[2],ncoord;
 char *dollar,actual_ind,subst;
 char temp[32] = "";
 char msg[80] = "";
 char tensor_name[32] = "";
 char tensor_1[32] = "";
 char tensor_2[32] = "";
 char t_indx[25] = "";
 char swap1,swap2,buffer[8192] = "";
 char seqbuf[8192] = "";
 char eq[2] = "";
 char swap[2] = "";
 char truInd[2] = "";

if (gr_debug) fprintf(stderr,"TRACE: %s %s %s\n",tensor,eq,expression);
 for (j=0;j<25;j++) 
    t_indx[j] = '\0';
 for (j=0;j<32;j++)
    {
     temp[j] = '\0'; 
     tensor_name[j] = '\0';
    }
 strcpy(tensor_1,tensor); strcpy(tensor_2,tensor);
 if (rank == 0) 
   {
    number_indices = 0;
    len = strlen(tensor_2);
    swap1 = tensor_2[len-1];
    swap2 = tensor_2[len-2];
    if (isupper(swap1))
      tensor_2[len-1] = toupper(swap2);
    else
      tensor_2[len-1] = tolower(swap2);
    if (isupper(swap2))
      tensor_2[len-2] = toupper(swap1);
    else
      tensor_2[len-2] = tolower(swap1);
   }
 else
   {
    mrank = is_index(tensor) + is_coord(tensor);
    dollar = (char *)index(tensor,'$');
    l = 0; n = 0;
    for (m=0; m<mrank; m++)
       {
	if (index( coord_list, *(dollar+m+1) ) )
	  {swap[l] = *(dollar+m+1); loc[l] = m; l++;}
	else
	  {truInd[n] = *(dollar+m+1); n++;}
       }
    dollar = (char *)index(tensor_2,'$');
    if (isupper(*(dollar+loc[0]+1)))
      *(dollar+loc[0]+1) = toupper(swap[1]);
    else
      *(dollar+loc[0]+1) = tolower(swap[1]);
    if (isupper(*(dollar+loc[1]+1)))
      *(dollar+loc[1]+1) = toupper(swap[0]);
    else
      *(dollar+loc[1]+1) = toupper(swap[0]);
    number_indices = parse_tensor(rank,tensor_2,tensor_name,t_indx);
    strcpy(temp,tensor_name);
    len = strlen(temp);
    temp[len] = '.';
    len++;
    temp[len] = '\0';
    strcat(temp,t_indx);    
   }
 strcpy(eq,"=");
 strcpy(seqbuf,seq); strcat(seqbuf,tensor_1); strcat(seqbuf," = ");
 strcat(seqbuf,expression); strcat(seqbuf," )");
 switch (number_indices)
       {
	case 0:
	       print_equation(tensor_2,eq,seqbuf);
	       break;
	case 1:
	       for (i=0;i<n_dim;i++)
		  {
		   strncpy(buffer,seqbuf,8192);
		   for (m=0;m<mrank;m++)
		      {
		       actual_ind = *(dollar+m+1);
		       if (index(index_list,actual_ind))
			 {
		          if (islower(actual_ind))
			    subst = lower_coord[i];
		          else
			    subst = upper_coord[i];
		          temp[len+m] = subst;
		          expression_replace(actual_ind,subst,buffer);
			 }
		      }
	           print_equation(temp,eq,buffer);
		  }
	       break;
	case 2:
	       for (i=0;i<n_dim;i++)
		  {
		   for (j=0;j<n_dim;j++)
		      {
		       strncpy(buffer,seqbuf,8192);
		       for (m=0;m<mrank;m++)
			  {
			   actual_ind = *(dollar+m+1);
			   if (index(index_list,actual_ind))
			     {
			      if (isupper(actual_ind))
			        {
			         if (actual_ind == toupper(truInd[0]))
				   subst = upper_coord[i];
			         else
				   subst = upper_coord[j];
			        }
			      else
			        {
			         if (actual_ind == tolower(truInd[0]))
				   subst = lower_coord[i];
			         else
				   subst = lower_coord[j];
			        }
			      temp[len+m] = subst;
			      expression_replace(actual_ind,subst,buffer);
			     }
			  }
		       print_equation(temp,eq,buffer);
		      }
		  }
	       break;
       }
}

/* Assign a tensor to an expression component by component */

void tensor_eq_expression(int rank, char *tensor, char *eq, char *expression)
{
 int  i,j,k,l,m,n,len,number_indices,mrank;
 char *dollar,actual_ind,subst;
 char temp[32] = "";
 char tend[32] = "";
 char tensor_name[32] = "";
 char t_indx[25] = "";
 char tindex[5] = "";
 char buffer[8192] = "";

 if (gr_debug) fprintf(stderr,"TRACE: %s %s %s\n",tensor,eq,expression);
 for (j=0;j<25;j++) 
    t_indx[j] = '\0';
 for (j=0;j<32;j++)
    {
     temp[j] = '\0'; 
     tensor_name[j] = '\0';
    }
 len = 0;
 if (rank > 0)
   {
    mrank = is_index(tensor) + is_coord(tensor);
    dollar = (char *)index(tensor,'$');
    strcpy(tend,(dollar+1));
    n = 0;
    for (m=0; m<mrank; m++)
       if (index(index_list, *(dollar+m+1)))
	 {tindex[n] = *(dollar+m+1); n++;}
    number_indices = parse_tensor(rank,tensor,tensor_name,t_indx);
    strcpy(temp,tensor_name);
    len = strlen(temp);
    temp[len] = '.';
    len++;
    temp[len] = '\0';
    strcat(temp,tend);
   }
 else
    number_indices = 0;
 switch (number_indices)
       {
	case 0:
	  print_equation(tensor,eq,expression);
	break;
	case 1:
	  for (i=0;i<n_dim;i++)
	     {
	      strncpy(buffer,expression,8192);
	      for (m=0;m<mrank;m++)
		 {
		  actual_ind = *(dollar+m+1);
		  if (index(index_list,actual_ind))
		    {
		     if (islower(actual_ind))
		       subst = lower_coord[i];
		     else
		       subst = upper_coord[i];
		     temp[len+m] = subst;
	             expression_replace(actual_ind,subst,buffer);
		    }
		 }
	      print_equation(temp,eq,buffer);
	     }
	  break;
	case 2:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  strncpy(buffer,expression,8192);
		  for (m=0;m<mrank;m++)
		     {
		      actual_ind = *(dollar+m+1);
		      if (index(index_list,actual_ind))
			{
		         if (isupper(actual_ind))
			   {
			    if (actual_ind == toupper(tindex[0]))
			      subst = upper_coord[i];
			    else
			      subst = upper_coord[j];
			   }
		         else
			   {
			    if (actual_ind == tolower(tindex[0]))
			      subst = lower_coord[i];
			    else
			      subst = lower_coord[j];
			   }
		         temp[len+m] = subst;
	                 expression_replace(actual_ind,subst,buffer);
			}
		     }
		  print_equation(temp,eq,buffer);
		 }
	     }
 	  break;
	case 3:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (k=0;k<n_dim;k++)
		     {
		      strncpy(buffer,expression,8192);
		      for (m=0;m<mrank;m++)
			 {
			  actual_ind = *(dollar+m+1);
			  if (index(index_list,actual_ind))
			    {
			     if (isupper(actual_ind))
			       {
			        if (actual_ind == toupper(tindex[0]))
			          subst = upper_coord[i];
			        else if (actual_ind == toupper(tindex[1]))
			          subst = upper_coord[j];
			        else
			          subst = upper_coord[k];
			       }
		             else
			       {
			        if (actual_ind == tolower(tindex[0]))
			          subst = lower_coord[i];
			        else if (actual_ind == tolower(tindex[1]))
			          subst = lower_coord[j];
			        else
			          subst = lower_coord[k];
			       }
		             temp[len+m] = subst;
	                     expression_replace(actual_ind,subst,buffer);
			    }
		         }
		      print_equation(temp,eq,buffer);
		     }
		 }
	     }
	  break;
	case 4:
	  for (i=0;i<n_dim;i++)
	     {
              for (j=0;j<n_dim;j++)
		 {
		  for (k=0;k<n_dim;k++)
		     {
		      for (l=0;l<n_dim;l++)
			 {
			  strncpy(buffer,expression,8192);
		          for (m=0;m<mrank;m++)
			     {
			      actual_ind = *(dollar+m+1);
			      if (index(index_list,actual_ind)) {
			      if (isupper(actual_ind))
			        {
			         if (actual_ind == toupper(tindex[0]))
			           subst = upper_coord[i];
			         else if (actual_ind == toupper(tindex[1]))
			           subst = upper_coord[j];
			         else if (actual_ind == toupper(tindex[2]))
			           subst = upper_coord[k];
			         else
			           subst = upper_coord[l];
			        }
		              else
			        {
			         if (actual_ind == tolower(tindex[0]))
			           subst = lower_coord[i];
			         else if (actual_ind == tolower(tindex[1]))
			           subst = lower_coord[j];
			         else if (actual_ind == tolower(tindex[2]))
			           subst = lower_coord[k];
			         else
			           subst = lower_coord[l];
			        }
		              temp[len+m] = subst;
	                      expression_replace(actual_ind,subst,buffer);
		             } }
	      		  print_equation(temp,eq,buffer);
			 }
		     }
		 }
	     }
	  break;
       }
}

/* Parse off the tensor into it's name, and indices */

int parse_tensor(int rank,char *tensor,char *tensor_name,char *t_indx)
{
 int i,j,number,true_numb;
 char tmp,up_tmp;
 char up_ten_ind[12] = "";

 for (j=0;j<12;j++) 
    {
     up_ten_ind[j] = '\0';
    }
 i = j = number = true_numb = 0;
 do {
     tensor_name[i] = tensor[i];
     i++;
    }
 while ( tensor[i] != '$');
 tensor_name[i] = '\0';
 do {
     i++;
     tmp = tensor[i];
     up_tmp = toupper(tmp);
     if ((isalpha(tmp)) && 
         (index(up_ten_ind,up_tmp) == '\0'))
       {
        t_indx[number] = tmp;
        up_ten_ind[number] = up_tmp;
        number++;
	if (index(upper_index,up_tmp)) true_numb++;
       }
    }
 while (tmp != '\0');
 t_indx[number] = '\0';
 return true_numb;
}

/* Get Tensor Name (prefix) */

void get_tensor_name(char *tensor, char *name)
{
 int i = 0;

 do {
     name[i] = tensor[i];
     i++;
    }
 while ( tensor[i] != '$');
 name[i] = '\0';
}

/* Get the rank of a tensor or an expression */

int get_rank(char *expression)
{
 int i;
 char free[25] = "";
 
 for (i=0;i<25;i++) free[i] = '\0';
 return (get_free_indices(expression,free));
}

/* Get the free indices from an expression and return rank */

int get_free_indices(char *expression, char *free)
{
 int i,j,ndumb,nindx;
 char dummy[25] = "";
 char indices[25] = "";

 for (i=0;i<25;i++)
    {
     dummy[i] = '\0';
     free[i] = '\0';
     indices[i] = '\0';
    }
 nindx = get_all_indices(expression,indices);
 ndumb = get_dummy_indices(expression,dummy);
 j = 0;
 for (i=0;i<nindx;i++)
    {
     if (index(dummy,indices[i]) == '\0')
       {
	free[j] = indices[i];
	j++;
       }
    }
 free[j] = '\0';
 return (strlen(free));
}

/* Get the dummy indices from an expression */

int get_dummy_indices(char *expression, char *dummy)
{
 int i,j,k,nindx;
 char opposite;
 char indices[25] = "";

 for (i=0;i<25;i++) indices[i] = '\0';
 nindx = get_all_indices(expression,indices);
 k = 0;
 for (i=0;i<nindx;i++)
    {
     if (isupper(indices[i]))
	opposite = tolower(indices[i]);
     else
	opposite = toupper(indices[i]);
     for (j=i+1;j<nindx;j++)
	{
	 if (indices[j] == opposite)
	   {
	    dummy[k] = indices[i];
	    k++;
	    dummy[k] = opposite;
	    k++;
	   }
	}
    }
 dummy[k] = '\0';
 return(strlen(dummy));
}

/* Get the all indices from a tensor */

int get_all_indices(char *expression, char *indices)
{
 int i;
 char *dollar,*tensor;
 char express[8192] = "";

 i = 0;
 strcpy(express,expression);
 tensor = strtok(express,"+-*/() ");
 while (tensor != NULL)
      {
       if ((dollar = (char *)index(tensor,'$')) != NULL)
	 {
	  dollar++;
	  while ( (*dollar != ' ')&&(*dollar != '\0') )
	       {
	        if ((index(indices,*dollar) == NULL) && 
		    (index(coord_list,*dollar) == NULL))
                  {
		   indices[i] = *dollar;
		   dollar++;
		   i++;
		  }
		else
		  dollar++;
               }
	 }
       tensor = strtok(NULL,"+-*/() ");
      }
 indices[i] = '\0';
 return (strlen(indices));
}

/* Compares two characters for the qsort routine */

int chcomp(char *ch1, char *ch2)
{
 if (*ch1 == *ch2)
   return 0;
 else if (*ch1 < *ch2)
   return -1;
 else
  return 1;
}

/* Compare all indices in two expression */

int cmp_all_indices(char *exp1, char *exp2)
{
 int  i,num_ind1,num_ind2,result;
 char free_ind1[25] = "";
 char free_ind2[25] = "";

 for (i=0;i<25;i++)
    {
     free_ind1[i] = '\0';
     free_ind2[i] = '\0';
    }
 num_ind1 = get_all_indices(exp1,free_ind1);
 num_ind2 = get_all_indices(exp2,free_ind2);
 qsort(free_ind1,num_ind1,sizeof(char),chcomp);
 qsort(free_ind2,num_ind2,sizeof(char),chcomp);
 if ( num_ind1 == num_ind2)
   {
    result = strcmp(free_ind1,free_ind2);
   }
 else
   {
    result = 1;
   }
 return result;
}

/* Compare the free indices in two tensors */

int cmp_free_indices(char *exp1, char *exp2)
{
 int  i,num_ind1,num_ind2,result;
 char free_ind1[25],free_ind2[25];

 for (i=0;i<25;i++)
    {
     free_ind1[i] = '\0';
     free_ind2[i] = '\0';
    }
 num_ind1 = get_free_indices(exp1,free_ind1);
 num_ind2 = get_free_indices(exp2,free_ind2);
 qsort(free_ind1,num_ind1,sizeof(char),chcomp);
 qsort(free_ind2,num_ind2,sizeof(char),chcomp);
 if ( num_ind1 == num_ind2)
   {
    result = strcmp(free_ind1,free_ind2);
   }
 else
   {
    result = 1;
   }
 return result;
}

/* Remove $'s from the expression and replace with _'s */

void remove_dollars(char* express)
{
 int i,length;

 length = strlen(express);
 for (i=0;i<length;i++)
    {
     if (express[i] == '$') express[i] = '.';
    }
}

/* Replace index character with coordinate character in expression */

void expression_replace(char indx, char coord, char *expression)
{
 int	i;
 char	tmp = '\0';

 i = 0;
 tmp = expression[i];
 while (tmp != '\0')
      {
       if (tmp == '$')
	 {
	  i++;
	  tmp = expression[i];
	  while (isalpha(tmp))
	       {
		if (tmp == indx) 
		  expression[i] = coord;
		i++;
		tmp = expression[i];
	       }
	 }
       else
	 {
          i++;
	  tmp = expression[i];
	 }
      }
}

/* negate an expression */

void negate_expression(char *expression)
{
 int	i;

 i = 0;
 if (expression[i] != '(')
   {
    while (expression[i] != '\0')
         {
          if (expression[i] == '+')
	    {
	     expression[i] = '-';
	     i++;
	    }
          else if (expression[i] == '-')
	    {
	     expression[i] = '+';
             i++;
	    }
          else
	    i++;
         }
   }
}

/* Print out an equation comtaining an expression in an intelligent fashion */

void print_equation(char *lhside, char *eq, char *rhside)
{
 int	i,j,n,indent,line_length,token_length,exp_length;
 char	bval,eval;
 char   spaces[80] = "";
 char   token[40] = "";

 indent = strlen(lhside) + 4;
 line_length = indent;
 for (n=0;n<indent;n++)
    spaces[n] = ' ';
 spaces[indent] = '\0';
 remove_dollars(rhside);
 printf(" %s %s ",lhside,eq);
 exp_length = strlen(rhside);
 for (i=0,j=1;j<=exp_length;j++)
    {
     if ((rhside[j] == '*')||(rhside[j]=='/')||
	 (rhside[j] == '+')||(rhside[j]=='-')||(rhside[j]=='\0'))
       {
	token_length = j - i;
        for (n=0;n<40;n++) token[n] = '\0';
	strncpy(token,&rhside[i],token_length);
	if ((line_length + token_length) > 79)
	  {
	   line_length = indent + token_length;
	   printf("\n%s",spaces);
	   no_lines_out++;
          }
	else
	  line_length += token_length;
	printf("%s",token);
	i = j;
       }
    }
 printf(";\n");
 no_lines_out++;
}

/* calculate covariant derivative */

void covariant_derivative(char *tensor, char* connect, char *deriv)
{
 int	i,j,len,nlen,clen,rank,orank;
 char	lo,hi,dummy,diff;
 char	otensor[32] = "";
 char   name[32] = "";
 char   cname[32] = "";
 char   ind[25] = "";
 char   oind[25] = "";
 char	temp[25] = "";
 char   tmp_name[32] = "";
 char   msg[80] = "";
 char   term[8192] = "";
 char   exp[8192] = "";

 if (strncmp(tensor,"_d_",3) != 0)
   {
    strcpy(msg,"First argument must be partial in covariant derivative");
    yyerror(msg); 
   }
 else
   {
    dummy = '\0';
    len = strlen(tensor);
    diff = tensor[len-1];
    rank = get_all_indices(tensor,ind);
    for (i=3,j=0;i<(len-1);i++)
       {
	if (tensor[i] != '$')
          {
           otensor[j] = tensor[i];
	   j++;
          }
	else
	  {
	   otensor[j-1] = '$';
	  }
       }
    otensor[j] = '\0';
    i = 0;
    while ((otensor[i] != '$')&&(otensor[i] != '\0'))
         {
	  name[i] = otensor[i];
          i++;
         }
    name[i] = '$';
    i++;
    name[i] = '\0';
    i = 0;
    while ((connect[i] != '$')&&(connect[i] != '\0'))
         {
	  cname[i] = connect[i];
          i++;
         }
    cname[i] = '$';
    i++;
    cname[i] = '\0';
    if (islower(connect[i])) 
      {
       strcpy(msg,"First index of connection coefficient must be upper");
       yyerror(msg);
      }
    orank = get_all_indices(otensor,oind);
    for (i=0;i<n_ind;i++)
       {
	lo = lower_index[i];
        hi = upper_index[i];
        if ((index(ind,lo) == '\0')&&(index(ind,hi) == '\0'))
          dummy = lo;
       }
    if (dummy == '\0')
      {
       strcpy(msg,"No dummy index available for covariant derivative");
       yyerror(msg); 
      }
    if (orank == 0)
      strcpy(deriv,tensor);
    else
      {
       strcpy(deriv,"( ");
       strcat(deriv,tensor);
      }
    nlen = strlen(name);
    clen = strlen(cname);
    for (i=0;i<orank;i++)
       {
	strcpy(temp,oind);
	if (isupper(oind[i]))
	  {
	   strcat(deriv," + ");
	   cname[clen] = toupper(oind[i]);
	   cname[clen+1] = tolower(dummy);
	   cname[clen+2] = diff;
	   cname[clen+3] = '\0';
	   temp[i] = toupper(dummy);
	   strcpy(tmp_name,name);
	   strcat(tmp_name,temp);
	   tensor_product(exp, cname, tmp_name);
	   strcat(deriv,exp);
	  }
	else
	  {
	   strcat(deriv," - ");
	   cname[clen] = toupper(dummy);
	   cname[clen+1] = tolower(oind[i]);
	   cname[clen+2] = diff;
	   cname[clen+3] = '\0';
	   temp[i] = tolower(dummy);
	   strcpy(tmp_name,name);
	   strcat(tmp_name,temp);
	   tensor_product(exp, cname, tmp_name);
	   strcat(deriv,exp);
	  }
       }
    if (orank > 0) strcat(deriv," )");
   }
 if (gr_debug)
   {
    fprintf(stderr,"TRACE: %s : %s ->\n",tensor,connect);
    fprintf(stderr,"       %s\n",deriv);
   }
}

/* calculate Lie derivative */

void lie_derivative(char *vector, char* tensor, char *deriv)
{
 int	i,j,vlen,tlen,dvlen,rank;
 char	lo,hi,dummy;
 char	dvector[32] = "";
 char   dtensor[32] = "";
 char   vname[32] = "";
 char   tname[32] = "";
 char   sd[2] = "";
 char   ind[25] = "";
 char	ttemp[32] = "";
 char   dvtemp[32] = "";
 char   msg[80] = "";
 char   exp[8192] = "";

 if (islower(vector[strlen(vector)-1])) 
   {
    strcpy(msg,"Lie Derivative direction must be upper index");
    yyerror(msg);
   }
 sd[1] = '\0';
 for (i=0; i<32; i++) dvector[i] = dtensor[i] = '\0';
 strcpy(dvector,"_d_");
 for (i=0,j=3;i<=strlen(vector);i++) {
    if (vector[i] != '$') {dvector[j] = vector[i]; j++;}
    else { dvector[j] = 'c'; j++;
           dvector[j] = '$'; j++; } }
 dvlen = strlen(dvector) - 1;
 strcpy(dtensor,"_d_");
 for (i=0,j=3;i<=strlen(tensor);i++) {
    if (tensor[i] != '$') {dtensor[j] = tensor[i]; j++;}
    else { dtensor[j] = 'c'; j++; dtensor[j] = '$'; j++;} }
 rank = get_all_indices(tensor,ind);
 for (i=0;i<n_ind;i++)
    {
     lo = lower_index[i];
     hi = upper_index[i];
     if ((index(ind,lo) == '\0')&&(index(ind,hi) == '\0')) dummy = lo;
    }
 if (dummy == '\0')
   {
    strcpy(msg,"No dummy index available for Lie derivative");
    yyerror(msg); 
   }
 sd[0] = dummy;
 strcat(dtensor,sd);
 i = 0;
 while ((vector[i] != '$')&&(vector[i] != '\0'))
      {
       vname[i] = vector[i];
       i++;
      }
 vname[i] = '$';
 i++;
 vname[i] = '\0';
 vlen = strlen(vname);
 i = 0;
 while ((tensor[i] != '$')&&(tensor[i] != '\0'))
      {
       tname[i] = tensor[i];
       i++;
      }
 tname[i] = '$';
 i++;
 tname[i] = '\0';
 tlen = strlen(tname);
 strcat(tname,ind);
 strcpy(deriv,"( ");
 vname[vlen] = toupper(dummy);
 vname[vlen+1] = '\0';
 tensor_product(exp,dtensor,vname);
 strcat(deriv,exp);
 strcpy(ttemp,tname);
 strcpy(dvtemp,dvector);
 if (rank > 0)
   {
    for (i=0; i<rank; i++)
       {
	strcpy(tname,ttemp);
	strcpy(dvector,dvtemp);
        if (isupper(ind[i]))
          {
	   strcat(deriv," - ");
	   tname[tlen+i] = toupper(dummy);
	   dvector[dvlen] = ind[i];
	   dvector[dvlen+1] = tolower(dummy);
	   dvector[dvlen+2] = '\0';
	   tensor_product(exp,tname,dvector);
	   strcat(deriv,exp);
          }
        else
          {
	   strcat(deriv," + ");
	   tname[tlen+i] = tolower(dummy);
	   dvector[dvlen] = toupper(dummy);
	   dvector[dvlen+1] = ind[i];
	   dvector[dvlen+2] = '\0';
	   tensor_product(exp,tname,dvector);
	   strcat(deriv,exp);
          }
       }
   }
 strcat(deriv," )");
 if (gr_debug)
   {
    fprintf(stderr,"TRACE: %s : %s ->\n",vector,tensor);
    fprintf(stderr,"       %s\n",deriv);
   }
}

/* expanded version of yyerror for grpp */

void yyerror(char *error_msg)
{
 fprintf(stderr," Syntax Error in Yacc Parser:\n");
 fprintf(stderr,"*** %s ***\n",error_msg);
 fprintf(stderr,"Error on input line: %d\n",no_lines_in);
}

/* write out the ANSI #line no and filename */

void write_line_no()
{
 if (line_flg)
   {
    printf("#line %d \"%s\" \n",no_lines_in,infile);
    no_lines_out++;
   }
}
