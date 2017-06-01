%{

/*
	File:		grpp_lex.l
	Description:	Lexical Analyzer for GRPP
	Author:		James Kent Blackburn
	Date:		Jan 2000

	Copyright (c) 1992-2000 
	By James Kent Blackburn
	All Rights Reserved
*/

#include <string.h>
#include <ctype.h>
#include "y.tab.h"

#define TRUE		1
#define FALSE		0

extern int	no_lines_in,no_lines_out,gr_debug;

int	n_dim,n_ind,
	grpp_okay = FALSE,
	grpp_block = FALSE,
	coord_flag = FALSE,
	index_flag = FALSE;

char	error_str[80] = "",
	coord_list[50] = "",
	upper_coord[25] = "",
	lower_coord[25] = "",
	index_list[50] = "",
	upper_index[25] = "",
	lower_index[25] = "",
	buffer[8192] = "";

void 	get_indices(char*),
	get_coordinates(char*),
	report_lex_error(char*),
	remove_derivative(char*, char*),
	expand_symmetry(char*);

int	is_vector(char*),
	is_index(char*),
	is_coord(char*);

%}

%a 3000
%o 4400

ws	[ \t]+
nl	\n
marker	("$$")
semi	[;]
equal	[+\-*/]?[=]
seq	("=+=")
aeq	("=-=")
block	[\{\}]
number	(([0-9]+)|([0-9]*\.[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9][0-9]*))
var	[*&]?[a-zA-Z_][a-zA-Z0-9_]*(\[[0-9]+\])*
func	[*&]?[a-zA-Z_][a-zA-Z0-9_]*[\(][a-zA-Z0-9_&\.\+\-\*\/(\), ]+[\)]
array	[*&]?[a-zA-Z][a-zA-Z0-9_]*\[\]
vector	(\{[a-zA-Z0-9_\.\+\-\*\/\(\) ]+[,][a-zA-Z0-9_\.\+\-\*\/\(\), ]+\})
typ	("rank1"|"rank2"|"rank3"|"rank4")[ \t]
r1	[*&]?[a-zA-Z][a-zA-Z0-9_]*\$[,;]?[a-zA-Z]
r2	[*&]?[a-zA-Z][a-zA-Z0-9_]*\$[,;]?[\[<]?[a-zA-Z][,;]?[a-zA-Z][\]>]?
r3	[*&]?[a-zA-Z][a-zA-Z0-9_]*\$[,;]?[\[<]?[a-zA-Z][,;]?[\[<]?[a-zA-Z][,;]?[\]>]?[a-zA-Z][\]>]?
r4	[*&]?[a-zA-Z][a-zA-Z0-9_]*\$[,;]?[\[<]?[a-zA-Z][,;]?[\[<]?[a-zA-Z][,;]?[\[<\]>]?[a-zA-Z][,;]?[\]>]?[a-zA-Z][\]>]?
comment	("/*"[^*/\n]*"*/")
coord	("coordinates$("[a-zA-Z0-9, ]+(");"|") ;"|")"))
indice	("indices$("[a-zA-Z, ]+(");"|") ;"|")"))
antisym	("totally_antisymmetric$")
covderv	("covariant_derivative$")
liederv	("Lie_derivative$")

%%

{ws}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		}
{semi}		{ if (!grpp_block) 
		    {
		     printf("%c",yytext[0]);
		    }
		  else
		    {
		     yylval.letter = yytext[0];
		     return SEMI;
		    }
		}
{equal}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,yytext);
		     return EQ;
		    }
		}
{seq}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,"( ");
		     return SYMEQ;
		    }
		}
{aeq}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,"-( ");
		     return SYMEQ;
		    }
		}
{block}		{ if (!grpp_block) 
		    {
		     printf("%c",yytext[0]);
		    }
		}
{nl}		{ no_lines_in++;
		  yylval.letter = yytext[0];
		  return NL; 
		}
{marker}	{ if (grpp_block)
		    {
		     grpp_block = FALSE;
		     strcpy(yylval.string,"/* Close GRPP Block */");
		    }
		  else
		    {
		     grpp_block = TRUE;
		     strcpy(yylval.string,"\n/* Open  GRPP Block */\n");
		     no_lines_out += 2;
		    }
		  return MARKER;
		}
{number}	{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,yytext);
		     return NUMBER;
		    }
		}
{var}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,yytext);
		     return VARIABLE;
		    }
		}
{func}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,yytext);
		     return FUNCTION;
		    }
		}
{array}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,yytext);
		     return ARRAY;
		    }
		}
{vector}	{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     if (is_vector(yytext))
		       {
			strcpy(yylval.string,yytext);
			return VECTOR;
		       }
		     else
		       {
			strcpy(error_str,"Vector definition Error");
			report_lex_error(error_str);
			return 0;			
		       }
		    }
		}
{typ}		{ if (!grpp_block) 
		    {
		     printf("%s",yytext);
		    }
		  else
		    {
		     strcpy(yylval.string,yytext);
		     return TYPE;
		    }
		}
{r1}		{ if ((is_index(yytext) == 1) && (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK1;
		    }
		  else if ((is_coord(yytext) == 1) && 
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     remove_dollars(buffer);
		     strcpy(yylval.string,buffer);
		     return VARIABLE;
		    }
		  else
		    {
		     strcpy(error_str,"Rank1 definition Error: ");
		     strcat(error_str,yytext);
		     report_lex_error(error_str);
		     return 0;
		    }
		}
{r2}		{ if ((is_index(yytext) == 2) && (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK2;
		    }
		  else if ((is_coord(yytext) == 2) && 
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     remove_dollars(buffer);
		     strcpy(yylval.string,buffer);
		     return VARIABLE;
		    }
		  else if ((is_index(yytext) == 1) && 
			   (is_coord(yytext) == 1) &&
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK1;
		    }
		  else
		    {
		     strcpy(error_str,"Rank2 definition Error: ");
		     strcat(error_str,yytext);
		     report_lex_error(error_str);
		     return 0;
		    }
		}
{r3}		{ if ((is_index(yytext) == 3) && (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK3;
		    }
		  else if ((is_coord(yytext) == 3) && 
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     remove_dollars(buffer);
		     strcpy(yylval.string,buffer);
		     return VARIABLE;
		    }
		  else if ((is_index(yytext) == 2) && 
			   (is_coord(yytext) == 1) &&
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK2;
		    }
		  else if ((is_index(yytext) == 1) && 
			   (is_coord(yytext) == 2) &&
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK1;
		    }
		  else
		    {
		     strcpy(error_str,"Rank3 definition Error: ");
		     strcat(error_str,yytext);
		     report_lex_error(error_str);
		     return 0;
		    }
		}
{r4}		{ if ((is_index(yytext) == 4) && (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK4;
		    }
		  else if ((is_coord(yytext) == 4) && 
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     remove_dollars(buffer);
		     strcpy(yylval.string,buffer);
		     return VARIABLE;
		    }
		  else if ((is_index(yytext) == 3) && 
			   (is_coord(yytext) == 1) &&
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK3;
		    }
		  else if ((is_index(yytext) == 2) && 
			   (is_coord(yytext) == 2) &&
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK2;
		    }
		  else if ((is_index(yytext) == 1) && 
			   (is_coord(yytext) == 3) &&
			   (grpp_block && grpp_okay))
		    {
		     remove_derivative(yytext, buffer);
		     strcpy(yylval.string,buffer);
		     return RANK1;
		    }
		  else
		    {
		     strcpy(error_str,"Rank4 definition Error: ");
		     strcat(error_str,yytext);
		     report_lex_error(error_str);
		     return 0;
		    }
		}
{comment}	{ strcpy(yylval.string,yytext);
		  return COMMENT; 
		}
{coord}		{ if (grpp_block)
		    {
		     get_coordinates(yytext);
		     coord_flag = TRUE;
		     if (index_flag) grpp_okay = TRUE;
		       return COORD;
		    }
		  else
		    {
		     strcpy(error_str,"coordinate$ Outside GRPP Block");
		     report_lex_error(error_str);
		     return 0;
		    }
		}
{indice}	{ if (grpp_block)
		    {
		     get_indices(yytext);
		     index_flag = TRUE;
		     if (coord_flag) grpp_okay = TRUE;
		       return INDICE;
		    }
		  else
		    {
		     strcpy(error_str,"indices$ Outside GRPP Block");
		     report_lex_error(error_str);
		     return 0;
		    }
		}
{antisym}	{ if (grpp_block) 
		    {
		     return ANTISYM;
		    }
		}
{covderv}	{ if (grpp_block) 
		    {
		     return COVDERIV;
		    }
		}
{liederv}	{ if (grpp_block) 
		    {
		     return LIEDERIV;
		    }
		}
.		{ if (!grpp_block) 
		    {
		     printf("%c",yytext[0]);
		    }
		  else
		    {
		     return yytext[0];
		    }
		}

%%

/* integer function to check that vector has right dimensions */

int is_vector(char *vector)
{
 int i,length,n_components;

 n_components = 1;
 length = strlen(vector);
 for (i=0;i<length;i++)
    if (vector[i] == ',') n_components++;
 if (n_components == n_dim)
   return 1;
 else
   return 0;
}

/* integer function which test for good index values in tensor */

int is_index(char *tensor)
{
 int   n;
 char *dollar,tmp,*check;

 n = 0;
 dollar = (char *)index(tensor,'$');
 do {
     tmp = *(++dollar);
     if (isalpha(tmp))
       {
        check = (char *)index(index_list,tmp);
        if (check != '\0') n++;
	else check = (char *)index(coord_list,tmp);
       }
     else
       check = (char *)index(",;[]<>",tmp);
     }
 while ((tmp != '\0') && (check != '\0'));
 return n;
}

/* integer function which test for good coordinate values in tensor */

int is_coord(char *tensor)
{
 int   n;
 char *dollar,tmp,*check;

 n = 0;
 dollar = (char *)index(tensor,'$');
 do {
     tmp = *(++dollar);
     if (isalpha(tmp))
       {
        check = (char *)index(coord_list,tmp);
        if (check != '\0') n++;
	else check = (char *)index(index_list,tmp);
       }
     else
       check = (char *)index(",;[]<>",tmp);
     }
 while ((tmp != '\0') && (check != '\0'));
 return n;
}

/* void function which reads coordinates out of coordinate statement */

void get_coordinates(char *coord_state)
{
 char *open_par,*close_par,tmp;
 int i;

 open_par  = (char *)index(coord_state,'(');
 close_par = (char *)rindex(coord_state,')');
 n_dim = 0;
 i = 0;
 do {
     tmp = *(++open_par);
     if ((isalnum(tmp)) && (index(coord_list,tmp) == '\0'))
       {
        coord_list[i] = toupper(tmp);
        upper_coord[n_dim] = toupper(tmp);
        i++;
	coord_list[i] = tolower(tmp);
        lower_coord[n_dim] = tolower(tmp);
	i++;
        n_dim++;
       }
    }
 while (open_par != close_par);
 coord_list[i] = '\0';
 upper_coord[n_dim] = '\0';
 lower_coord[n_dim] = '\0';
}

/* void function which reads indices out of index statement */

void get_indices(char *index_state)
{
 char *open_par,*close_par,tmp;
 int i;

 open_par  = (char *)index(index_state,'(');
 close_par = (char *)rindex(index_state,')');
 n_ind = 0;
 i = 0;
 do {
     tmp = *(++open_par);
     if ((isalnum(tmp)) && (index(index_list,tmp) == '\0'))
       {
        index_list[i] = toupper(tmp);
        upper_index[n_ind] = toupper(tmp);
        i++;
	index_list[i] = tolower(tmp);
        lower_index[n_ind] = tolower(tmp);
	i++;
        n_ind++;
       }
    }
 while (open_par != close_par);
 index_list[i] = '\0';
 upper_index[n_ind] = '\0';
 lower_index[n_ind] = '\0';
}

/* remove derivative syntax (,;) from tensor string */

void remove_derivative(char *tensor, char *new_tensor)
{
 int	i,j,k,partial,covariant,length,order,ptr;
 char	clean[256] = "";

 length = strlen(tensor);
 j = ptr = partial = covariant = order = 0;
 for (i=0;i<length;i++)
    {
     if (tensor[i] == ',')
       partial = 1;
     else if (tensor[i] == ';')
       covariant = 1;
     else if (isalpha(tensor[i]))
       {
	clean[j] = tensor[i];
	j++;
	if (partial||covariant) order++;
       }
     else
       {
        if ((i==0) && (j == 0) && (tensor[i] == '*'))
          {
	   ptr = 1;
          }
        else if (tensor[i] == '$')
          { 
           clean[j] = '_'; j++;
	   for (k=(i+1);k<length;k++)
              {
	       if (isalpha(tensor[k]))
		 {
	          if (isupper(tensor[k])) 
		    clean[j] = 'C';
	          else
		    clean[j] = 'c';
	          j++;
		 }
              }
	   if (ptr)
	     {
	      clean[j] = ')';
	      j++;
	     }
	  }
	if (tensor[i] != '*') {
          clean[j] = tensor[i];
	  j++; }
       }
    }
 clean[j] = '\0';
 if (ptr) strcpy(new_tensor,"(*");
 else new_tensor[0] = '\0';
 switch (order)
       {
	case 0:
	  strcat(new_tensor,clean);
	break;
	case 1:
	  if (partial)
	    strcat(new_tensor,"_d_");
	  else if (covariant)
	    strcat(new_tensor,"_D_");
	  strcat(new_tensor,clean);
	break;
	case 2:
	  if (partial)
	    strcat(new_tensor,"_d2_");
	  else if (covariant)
	    strcat(new_tensor,"_D2_");
	  strcat(new_tensor,clean);
	break;
	case 3:
	  if (partial)
	    strcat(new_tensor,"_d3_");
	  else if (covariant)
	    strcat(new_tensor,"_D3_");
	  strcat(new_tensor,clean);
	break;
	case 4:
	  if (partial)
	    strcat(new_tensor,"_d4_");
	  else if (covariant)
	    strcat(new_tensor,"_D4_");
	  strcat(new_tensor,clean);
	break;
       }
 if (gr_debug)
   {
    fprintf(stderr,"TRACE: %s ->\n",tensor);
    fprintf(stderr,"       %s\n",new_tensor);
   }
 expand_symmetry(new_tensor);
}

/* expand any anti/symmetries */

void expand_symmetry(char* tensor)
{
 int	i,j,k,l,n,value,error,elements,sign,total,factor,nlist[5];
 char	*opn,list[5];
 char   ssign[5] = "";
 char   sfactor[5] = "";
 char	msg[80] = "";
 char   name[32] = "";
 char   expression[8192] = "";

 for (i=0;i<5;i++) list[i] = '\0';
 error = elements = sign = n = 0;
 if ((opn = (char *)index(tensor,'<')) != '\0')
   {
    if (index(tensor,'>') != '\0')
      {
       while (tensor[n] != '<')
	    {
	     name[n] = tensor[n];
	     n++;
	    }
       name[n] = '\0';
       sign = +1;
       while ((*opn != '>')&&(*opn != '\0'))
            {
	     if (isalpha(*opn))
	       {
		list[elements] = *opn;
		elements++;
	       }
	     opn++;
            }
       list[elements] = '\0';
       if ((elements < 2)||(elements > 4)) error = 1;
      }
    else
      {
       error = 1;
       strcpy(msg,"Error in symmetric tensor operator");
      }
   }
 else if ((opn = (char *)index(tensor,'[')) != '\0')
   {
    if (index(tensor,']') != '\0')
      {
       while (tensor[n] != '[')
	    {
	     name[n] = tensor[n];
	     n++;
	    }
       name[n] = '\0';
       sign = -1;
       while ((*opn != ']')&&(*opn != '\0'))
            {
	     if (isalpha(*opn))
	       {
		list[elements] = *opn;
		elements++;
	       }
	     opn++;
            }
       list[elements] = '\0';
       if ((elements < 2)||(elements > 4)) error = 1;
      }
    else
      {
       error = 1;
       strcpy(msg,"Error in anti-symmetric tensor operator");
      }
   }
 else if (index(tensor,'>') != '\0')
   {
    error = 1;
    strcpy(msg,"Error in symmetric tensor operator");
   }
 else if (index(tensor,']') != '\0')
   {
    error = 1;
    strcpy(msg,"Error in anti-symmetric tensor operator");
   }
 else return;
 if (error) report_lex_error(msg);
 if (sign != 0)
   {
    total = 0;
    strcpy(expression,"( ( ");
    switch (elements)
	  {
	   case 2:
	     factor = 2;
	     strcpy(sfactor,"2");
	     for (i=0;i<elements;i++)
		for (j=0;j<elements;j++)
		   {
		    if ( !(i==j) )
		      {
		       if (sign == -1)
			 {
		          nlist[0] = i; nlist[1] = j;
		          value = permutation_sort(elements,nlist);
		          if (value == 1)
			    strcpy(ssign," + ");
		          else
			    strcpy(ssign," - ");
			 }
		       else
			 strcpy(ssign," + ");
		       name[n] = list[i];
		       name[n+1] = list[j];
		       name[n+2] = '\0';
		       if ( (total > 0)&&(total < factor) )
			 strcat(expression,ssign);
		       strcat(expression,name);
		       total++;
		      }
		   }
	   break;
	   case 3:
	     factor = 6;
	     strcpy(sfactor,"6");
	     for (i=0;i<elements;i++)
		for (j=0;j<elements;j++)
		   for (k=0;k<elements;k++)
		      {
		       if ( !((i==j)||(i==k)||(j==k)) )
		         {
		          if (sign == -1)
			    {
		             nlist[0] = i; nlist[1] = j; nlist[2] = k;
		             value = permutation_sort(elements,nlist);
		             if (value == 1)
			       strcpy(ssign," + ");
		             else
			       strcpy(ssign," - ");
			    }
		          else
			    strcpy(ssign," + ");
		          name[n] = list[i];
		          name[n+1] = list[j];
		          name[n+2] = list[k];
		          name[n+3] = '\0';
		          if ( (total > 0)&&(total < factor) )
			    strcat(expression,ssign);
		          strcat(expression,name);
		          total++;
		         }
		      }
	   break;
	   case 4:
	     factor = 24;
	     strcpy(sfactor,"24");
	     for (i=0;i<elements;i++)
		for (j=0;j<elements;j++)
		   for (k=0;k<elements;k++)
		      for (l=0;l<elements;l++)
		         {
		          if ( !((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) )
		            {
		             if (sign == -1)
			       {
		                nlist[0] = i; nlist[1] = j;
			        nlist[2] = k; nlist[3] = l;
		                value = permutation_sort(elements,nlist);
		                if (value == 1)
			          strcpy(ssign," + ");
		                else
			          strcpy(ssign," - ");
			       }
		             else
			       strcpy(ssign," + ");
		             name[n] = list[i];
		             name[n+1] = list[j];
		             name[n+2] = list[k];
		             name[n+3] = list[l];
		             name[n+4] = '\0';
		             if ( (total > 0)&&(total < factor) )
			       strcat(expression,ssign);
		             strcat(expression,name);
		             total++;
		            }
		         }
	   break;
	  }
    strcat(expression," )");
    strcat(expression," / ");
    strcat(expression,sfactor);
    strcat(expression," )");
    if ( (index(expression,'<') != '\0') || (index(expression,'<') != '\0') ||
         (index(expression,'[') != '\0') || (index(expression,']') != '\0') )
      {
       strcpy(msg,"Error in anti/symmetry operator: two occurrences?");
       report_lex_error(msg);
      }
    else
      {
       if (gr_debug)
	 {
	  fprintf(stderr,"TRACE: %s ->\n",tensor);
	  fprintf(stderr,"       %s\n",expression);
	 }
       strcpy(tensor,expression);
      }
   }
}

/* report a lexical error to the user */

void report_lex_error(char *error_str)
{
 fprintf(stderr,"Error in Lexical Analyzer:\n");
 fprintf(stderr,"*** %s ***\n",error_str);
 fprintf(stderr,"Error on input line: %d\n",no_lines_in);
}
