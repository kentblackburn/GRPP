/*
	File:		grpp.c
	Description:	General Relativity C preprocessor
	Author:		James Kent Blackburn
	Date:		Jan 2000

	Copyright (c) 1992-2000 
	By James Kent Blackburn
	All Rights Reserved
*/

char *program_name = "grpp",
     *usage = "USAGE:\n %s [-trace] [-lines] [-o outfile] infile\n",
     *infile, *outfile, *hdrfile, hdrf[80];
int   no_lines_in  = 1,
      no_lines_out = 0; 
int   syntax_error = 0,
      gr_debug = 0,
      line_flg = 0;

#include <stdio.h>
#include <string.h>
#include <sysexits.h>

#define DEFAULT_OUTFILE "grpp_out.c"
#define DEFAULT_HDRFILE "grpp_out.h"

main (argc,argv)
 int argc;
 char **argv;

{
 int i=0,j;
 FILE *in,
      *out;
 void start(),finish();

 start();
 program_name = argv[0];
 infile = NULL;
 if (argc == 1)
   {
    fprintf(stderr,usage,program_name);
    exit(EX_USAGE);
   }
 else
   {
    outfile = DEFAULT_OUTFILE;
    hdrfile = DEFAULT_HDRFILE;
    for (j=1; j<argc; j++)
       {
	if (!strcmp(argv[j],"-trace"))
	  gr_debug = 1;
	else if (!strcmp(argv[j],"-lines"))
	  line_flg = 1;
	else if (!strcmp(argv[j],"-o"))
	  {
	   j++;
	   outfile = argv[j];
	   do { hdrf[i] = outfile[i];
	        i++;
	      } while ( (outfile[i] != '.') && (outfile[i] != '\0'));
	   hdrf[i] = '\0';
	   strcat(hdrf,"_grpp.h");
	   hdrfile = &hdrf[0];
	  }
	else
	  infile = argv[j];
       }
    if (gr_debug) fprintf(stderr,"Trace enabled\n");
    if (line_flg) fprintf(stderr,"#Line enabled\n");
    if (infile == NULL)
      {
       fprintf(stderr,"No infile specified\n");
       fprintf(stderr,usage,program_name);
       exit(EX_USAGE);
      }
    in = freopen(infile,"r",stdin);
    if (in == NULL)
      {
       fprintf(stderr,"Cannot find input: %s\n",infile);
       exit(EX_NOINPUT);
      }
    out = freopen(outfile,"w",stdout);
    if (out == NULL)
      {
       fprintf(stderr,"Cannot create output: %s\n",outfile);
       exit(EX_CANTCREAT);
      }
    if (line_flg)
      {
       printf("#line %d \"%s\" \n",no_lines_in,infile);
       no_lines_out++;
      }
    yyparse();
    fclose(out);
    fclose(in);
    if (syntax_error)
      {
       fprintf(stderr,"Syntax error occurred \n");
       unlink(outfile);
       exit(EX_DATAERR);
      }
    else
      {
       finish(infile,outfile);
       exit(0);
      }
   }
}

void start()
{
 fprintf(stderr,"\nGeneral Relativity PreProcessor \n");
 fprintf(stderr,"         Version 2.1 \n\n");
 fprintf(stderr,"     Copyright (c) 1992-2000 \n");
 fprintf(stderr,"     by James Kent Blackburn \n");
 fprintf(stderr,"     All Rights Reserved \n\n");
}

void finish(infile,outfile)
 char *infile, *outfile;
{
 fprintf(stderr,"Results: \n");
 fprintf(stderr,"     %s: %d total lines \n",infile, no_lines_in );
 fprintf(stderr,"     %s: %d total lines \n\n",outfile,no_lines_out);
}


int yywrap()
{
 return ( 1 );
}


