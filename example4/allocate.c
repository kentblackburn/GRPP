#include <stdlib.h>
#include <stdio.h>

void nrerror(error_text)
char error_text[];
{
 void exit();

 fprintf(stderr,"Numerical Recipes Run-Time Error:\n");
 fprintf(stderr,"%s\n",error_text);
 fprintf(stderr,"Now Exiting to system\n");
 exit(1);
}

double *vector(nl,nh)
int nl,nh;
{
 double *v;
 
 v = (double *) malloc((unsigned) (nh-nl+1) * sizeof(double));
 if (!v) nrerror("allocation failure in vector()");
 return v-nl;
}

double **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
 int i;
 double **m;

 m = (double **) malloc((unsigned) (nrh-nrl+1) * sizeof(double*));
 if (!m) nrerror("row allocation failure in matrix(r,c)");
 m -= nrl;
 for (i=nrl;i<=nrh;i++)
    {
     m[i] = (double *) malloc((unsigned) (nch-ncl+1) * sizeof(double));
     if (!m[i]) nrerror("column allocation failure in matrix(r,c)");
     m[i] -= ncl;
    }
 return m;
}

void free_vector(v,nl,nh)
double *v;
int nl,nh;
{
 free((char *) (v+nl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
 int i;

 for (i=nrh;i>=nrl;i--) free((char *) (m[i] + ncl));
 free((char *) (m + nrl));
}
