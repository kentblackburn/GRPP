%{

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

%}

%union	{
	 char	letter;
	 char	string[8192];
	}

%token	<letter>	NL SEMI
%token	<string>	EQ SYMEQ MARKER NUMBER VARIABLE FUNCTION
%token	<string>	TYPE RANK1 RANK2 RANK3 RANK4 COORD INDICE
%token	<string>	ARRAY VECTOR ANTISYM COMMENT COVDERIV LIEDERIV

%type	<string>	expression

%right	EQ SYMEQ
%left	'+' '-'
%right	'*' '/'
%left	UMINUS

%%

lines:		  /* nothing */
		| lines line
		;
		
line:		NL
			{no_lines_out++; printf("\n");
			 write_line_no();}
		| MARKER
			{printf("%s",$1);}
		| COMMENT
			{printf("%s",$1);}
		| COORD
			{printf("/* COORDINATES: %s */\n",coord_list);
			 printf("/* DIMENSION: %d */\n",n_dim);
			 printf("\n#include \"%s\" \n",hdrfile);
			 no_lines_out += 4;
			 tensor_structures();}
		| INDICE
			{printf("/* INDICES: %s */\n",index_list);
			 printf("/* TOTAL: %d */\n",n_ind);
			 no_lines_out += 2;
			 tensor_structures();}
		| declare
		| assign
		;

declare:	 TYPE RANK1 SEMI
			{tensor_declare(1,$1,$2);}
		| TYPE RANK2 SEMI
			{tensor_declare(2,$1,$2);}
		| TYPE RANK3 SEMI
			{tensor_declare(3,$1,$2);}
		| TYPE RANK4 SEMI
			{tensor_declare(4,$1,$2);}
		| TYPE RANK1 ',' RANK1 SEMI
			{tensor_declare(1,$1,$2);
			 tensor_declare(1,$1,$4);}
		| TYPE RANK2 ',' RANK2 SEMI
			{tensor_declare(2,$1,$2);
			 tensor_declare(2,$1,$4);}
		| TYPE RANK3 ',' RANK3 SEMI
			{tensor_declare(3,$1,$2);
			 tensor_declare(3,$1,$4);}
		| TYPE RANK4 ',' RANK4 SEMI
			{tensor_declare(4,$1,$2);
			 tensor_declare(4,$1,$4);}
		| TYPE RANK1 ',' RANK1 ',' RANK1 SEMI
			{tensor_declare(1,$1,$2);
			 tensor_declare(1,$1,$4);
			 tensor_declare(1,$1,$6);}
		| TYPE RANK2 ',' RANK2 ',' RANK2 SEMI
			{tensor_declare(2,$1,$2);
			 tensor_declare(2,$1,$4);
			 tensor_declare(2,$1,$6);}
		| TYPE RANK3 ',' RANK3 ',' RANK3 SEMI
			{tensor_declare(3,$1,$2);
			 tensor_declare(3,$1,$4);
			 tensor_declare(3,$1,$6);}
		| TYPE RANK4 ',' RANK4 ',' RANK4 SEMI
			{tensor_declare(4,$1,$2);
			 tensor_declare(4,$1,$4);
			 tensor_declare(4,$1,$6);}
		| TYPE RANK1 ',' RANK1 ',' RANK1 ',' RANK1 SEMI
			{tensor_declare(1,$1,$2);
			 tensor_declare(1,$1,$4);
			 tensor_declare(1,$1,$6);
			 tensor_declare(1,$1,$8);}
		| TYPE RANK2 ',' RANK2 ',' RANK2 ',' RANK2 SEMI
			{tensor_declare(2,$1,$2);
			 tensor_declare(2,$1,$4);
			 tensor_declare(2,$1,$6);
			 tensor_declare(2,$1,$8);}
		| TYPE RANK3 ',' RANK3 ',' RANK3 ',' RANK3 SEMI
			{tensor_declare(3,$1,$2);
			 tensor_declare(3,$1,$4);
			 tensor_declare(3,$1,$6);
			 tensor_declare(3,$1,$8);}
		| TYPE RANK4 ',' RANK4 ',' RANK4 ',' RANK4 SEMI
			{tensor_declare(4,$1,$2);
			 tensor_declare(4,$1,$4);
			 tensor_declare(4,$1,$6);
			 tensor_declare(4,$1,$8);}
		| TYPE RANK1 ',' RANK1 ',' RANK1 ',' RANK1 ',' RANK1 SEMI
			{tensor_declare(1,$1,$2);
			 tensor_declare(1,$1,$4);
			 tensor_declare(1,$1,$6);
			 tensor_declare(1,$1,$8);
			 tensor_declare(1,$1,$10);}
		| TYPE RANK2 ',' RANK2 ',' RANK2 ',' RANK2 ',' RANK2 SEMI
			{tensor_declare(2,$1,$2);
			 tensor_declare(2,$1,$4);
			 tensor_declare(2,$1,$6);
			 tensor_declare(2,$1,$8);
			 tensor_declare(2,$1,$10);}
		| TYPE RANK3 ',' RANK3 ',' RANK3 ',' RANK3 ',' RANK3 SEMI
			{tensor_declare(3,$1,$2);
			 tensor_declare(3,$1,$4);
			 tensor_declare(3,$1,$6);
			 tensor_declare(3,$1,$8);
			 tensor_declare(3,$1,$10);}
		| TYPE RANK4 ',' RANK4 ',' RANK4 ',' RANK4 ',' RANK4 SEMI
			{tensor_declare(4,$1,$2);
			 tensor_declare(4,$1,$4);
			 tensor_declare(4,$1,$6);
			 tensor_declare(4,$1,$8);
			 tensor_declare(4,$1,$10);}
		| TYPE RANK1 ',' RANK1 ',' RANK1 ',' RANK1 ',' RANK1 ',' RANK1 SEMI
			{tensor_declare(1,$1,$2);
			 tensor_declare(1,$1,$4);
			 tensor_declare(1,$1,$6);
			 tensor_declare(1,$1,$8);
			 tensor_declare(1,$1,$10);
			 tensor_declare(1,$1,$12);}
		| TYPE RANK2 ',' RANK2 ',' RANK2 ',' RANK2 ',' RANK2 ',' RANK2 SEMI
			{tensor_declare(2,$1,$2);
			 tensor_declare(2,$1,$4);
			 tensor_declare(2,$1,$6);
			 tensor_declare(2,$1,$8);
			 tensor_declare(2,$1,$10);
			 tensor_declare(2,$1,$12);}
		| TYPE RANK3 ',' RANK3 ',' RANK3 ',' RANK3 ',' RANK3 ',' RANK3 SEMI
			{tensor_declare(3,$1,$2);
			 tensor_declare(3,$1,$4);
			 tensor_declare(3,$1,$6);
			 tensor_declare(3,$1,$8);
			 tensor_declare(3,$1,$10);
			 tensor_declare(3,$1,$12);}
		| TYPE RANK4 ',' RANK4 ',' RANK4 ',' RANK4 ',' RANK4 ',' RANK4 SEMI
			{tensor_declare(4,$1,$2);
			 tensor_declare(4,$1,$4);
			 tensor_declare(4,$1,$6);
			 tensor_declare(4,$1,$8);
			 tensor_declare(4,$1,$10);
			 tensor_declare(4,$1,$12);}
		;

assign:		  ANTISYM '(' VARIABLE ')' SEMI
			{total_antisymmetric_tensor($3);}
		| FUNCTION SEMI
			{printf(" %s;\n",$1); no_lines_out++;}
		| RANK1 EQ VECTOR SEMI
			{vector_equation($1,$2,$3);}
		| VECTOR EQ RANK1 SEMI
			{equation_vector($1,$2,$3);}
		| RANK1 EQ ARRAY SEMI
			{tensor_eq_array(1,$1,$2,$3);}
		| RANK2 EQ ARRAY SEMI
			{tensor_eq_array(2,$1,$2,$3);}
		| RANK3 EQ ARRAY SEMI
			{tensor_eq_array(3,$1,$2,$3);}
		| RANK4 EQ ARRAY SEMI
			{tensor_eq_array(4,$1,$2,$3);}
		| ARRAY EQ expression SEMI
			{array_eq_expression($1,$2,$3);}
		| VARIABLE EQ expression SEMI
			{tensor_eq_expression(0,$1,$2,$3);}
		| RANK1 EQ expression SEMI
			{tensor_eq_expression(1,$1,$2,$3);}
		| RANK2 EQ expression SEMI
			{tensor_eq_expression(2,$1,$2,$3);}
		| RANK3 EQ expression SEMI
			{tensor_eq_expression(3,$1,$2,$3);}
		| RANK4 EQ expression SEMI
			{tensor_eq_expression(4,$1,$2,$3);}
		| VARIABLE SYMEQ expression SEMI
			{symten_eq_expression(0,$1,$2,$3);}
		| RANK1 SYMEQ expression SEMI
			{symten_eq_expression(1,$1,$2,$3);}
		| RANK2 SYMEQ expression SEMI
			{symten_eq_expression(2,$1,$2,$3);}
		;

expression:	  NUMBER
			{
			 strcpy($$,$1);
			}
		| VARIABLE
			{
			 strcpy($$,$1);
			}
		| FUNCTION
			{
			 strcpy($$,$1);
			}
		| RANK1
			{
			 expand_dummy_indices(1,$1);
			 strcpy($$,$1);
			}
		| RANK2
			{
			 expand_dummy_indices(2,$1);
			 strcpy($$,$1);
			}
		| RANK3
			{
			 expand_dummy_indices(3,$1);
			 strcpy($$,$1);
			}
		| RANK4
			{
			 expand_dummy_indices(4,$1);
			 strcpy($$,$1);
			}
		| COVDERIV '(' RANK1 ':' RANK3 ')' 
			{ covariant_derivative($3,$5,$$); }
		| COVDERIV '(' RANK2 ':' RANK3 ')' 
			{ covariant_derivative($3,$5,$$); }
		| COVDERIV '(' RANK3 ':' RANK3 ')' 
			{ covariant_derivative($3,$5,$$); }
		| COVDERIV '(' RANK4 ':' RANK3 ')' 
			{ covariant_derivative($3,$5,$$); }
		| LIEDERIV '(' RANK1 ':' RANK1 ')' 
			{ lie_derivative($3,$5,$$); }
		| LIEDERIV '(' RANK1 ':' RANK2 ')' 
			{ lie_derivative($3,$5,$$); }
		| LIEDERIV '(' RANK1 ':' RANK3 ')' 
			{ lie_derivative($3,$5,$$); }
		| LIEDERIV '(' RANK1 ':' RANK4 ')' 
			{ lie_derivative($3,$5,$$); }
		| expression '+' expression
			{
			 if (!sum_expressions($$, $1, '+', $3))
			   {
			    strcpy(msg,"Tensor Summation Error");
			    strcat(msg,$1);
			    strcat(msg," + ");
			    strcat(msg,$3);
			    yyerror(msg);
			   }
			}
		| expression '-' expression
			{
			 if (!sum_expressions($$, $1, '-', $3))
			   {
			    strcpy(msg,"Tensor Subtraction Error");
			    strcat(msg,$1);
			    strcat(msg," - ");
			    strcat(msg,$3);
			    yyerror(msg);
			   }
			}
		| expression '*' expression
			{
			 if (!tensor_product($$,$1,$3))
			   {
			    strcpy(msg,"Tensor Product Error");
			    strcat(msg,$1);
			    strcat(msg," * ");
			    strcat(msg,$3);
			    yyerror(msg);
			   }
			}
		| expression '/' expression
			{
			    strcpy($$,$1);
			    strcat($$," / ");
			    strcat($$,$3);
			}
		| '-' expression %prec UMINUS
			{
			 strcpy($$,"-");
			 strcat($$,$2);
			}
		| '(' expression ')'
			{
			 if (($2[0] == '(')&&($2[strlen($2)-1] == ')'))
			   strcpy($$,$2);
			 else
			   {
			    strcpy($$,"( ");
			    strcat($$,$2);
			    strcat($$," )");
			   }
			}
		;

%%

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
