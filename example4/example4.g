/* 
   	File:		example4.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Determine curvature for Kerr black hole

	Copyright (c) 1993-1998 
	By James Kent Blackburn
	All Rights Reserved
*/

$$ /* set up coordinates and indices used by GRPP */
coordinates$(t,r,a,p);
indices$(i,j,k,l,m);
$$

#include <stdio.h>
#include <math.h>

#define MASS	1.00	/* Mass of rotating black hole */
#define ANGM	0.99	/* Angular momentum per mass of rotating black hole */
#define EPS	1.0e-4	/* Epsilon, a small number */
#define PI	3.1415926535897932384626433

double **array4x4; /* global 4x4 array for inverse metric calls */
void inverse(double**,int);
double **matrix(int,int,int,int);

void kerr_metric(radius, theta, phi, g_cc)
double radius, theta, phi;
/* last argument was Tensor */
$$ rank2 *g$ij ; $$
{
 double r2, A2, cos2a, sin2a, sigma, delta;
 double gtt, grr, gaa, gpp, gtp;

/* build up common terms in kerr metric */
 r2 = radius * radius; A2 = ANGM * ANGM; 
 cos2a = cos( theta ); cos2a *= cos2a;
 sin2a = sin( theta ); sin2a *= sin2a;
 sigma = r2 + A2 * cos2a; delta = r2 - 2 * MASS * radius + A2;

/* do scalar algebra here, C is better at that! */
 gtt = - (1.0 - 2 * MASS * radius / sigma);
 grr = sigma / delta;
 gaa = sigma;
 gpp = (r2 + A2 + 2 * MASS * radius * A2 * sin2a / sigma) * sin2a;
 gtp = - 2 * ANGM * MASS * radius * sin2a / sigma;

$$ /* assign all components of kerr metric tensor */
 *g$tt = gtt;
 *g$rr = grr;
 *g$aa = gaa;
 *g$pp = gpp;
 *g$tp =+= gtp;
 *g$tr =+= 0.0;
 *g$ta =+= 0.0;
 *g$ra =+= 0.0;
 *g$rp =+= 0.0;
 *g$ap =+= 0.0;
$$
}

void kerr_christoffel(radius, theta, phi, Christoffel_Ccc)
double radius, theta, phi;
$$ rank3 *Christoffel$Ijk; $$
{
$$
 rank2 g$ij, g$IJ;
 rank2 g_r_plus_eps$ij, g_r_minus_eps$ij ;
 rank2 g_a_plus_eps$ij, g_a_minus_eps$ij ;
 rank3 g$ij,k;
$$
 void kerr_metric(), inverse();

 kerr_metric(radius, theta, phi, &g_cc);
$$ array4x4[] = g$ij; $$
 inverse(array4x4,4);
$$ g$IJ = array4x4[]; $$
 radius += EPS;
 kerr_metric(radius, theta, phi, &g_r_plus_eps_cc);
 radius -= 2*EPS;
 kerr_metric(radius, theta, phi, &g_r_minus_eps_cc);
 radius += EPS;
 theta  += EPS;
 kerr_metric(radius, theta, phi, &g_a_plus_eps_cc);
 theta  -= 2*EPS;
 kerr_metric(radius, theta, phi, &g_a_minus_eps_cc);
 theta  += EPS;
$$
 g$ij,t = 0.0;
 g$ij,r = ( g_r_plus_eps$ij - g_r_minus_eps$ij ) / (2 * EPS);
 g$ij,a = ( g_a_plus_eps$ij - g_a_minus_eps$ij ) / (2 * EPS);
 g$ij,p = 0.0;
 *Christoffel$Ijk = 0.5 * g$IL * ( g$lj,k + g$lk,j - g$kj,l );
$$
}

double kerr_curvature(radius, theta, phi)
double radius, theta, phi;
{
 double R;
$$
 rank2 g$ij, g$IJ, R$ij;
 rank3 Chr$Ijk;
 rank3 Chr_r_plus_eps$Ijk, Chr_r_minus_eps$Ijk;
 rank3 Chr_a_plus_eps$Ijk, Chr_a_minus_eps$Ijk;
 rank4 R$Ijkl, Chr$Ijk,l;
$$
 void kerr_metric(), kerr_christoffel(), inverse();

 kerr_metric(radius, theta, phi, &g_cc);
$$ array4x4[] = g$ij; $$
 inverse(array4x4,4);
$$ g$IJ = array4x4[]; $$
 kerr_christoffel(radius, theta, phi, &Chr_Ccc);
 radius += EPS;
 kerr_christoffel(radius, theta, phi, &Chr_r_plus_eps_Ccc);
 radius -= 2*EPS;
 kerr_christoffel(radius, theta, phi, &Chr_r_minus_eps_Ccc);
 radius += EPS;
 theta  += EPS;
 kerr_christoffel(radius, theta, phi, &Chr_a_plus_eps_Ccc);
 theta  -= 2*EPS;
 kerr_christoffel(radius, theta, phi, &Chr_a_minus_eps_Ccc);
 theta  += EPS;
$$
 Chr$Ijk,t = 0.0;
 Chr$Ijk,r = ( Chr_r_plus_eps$Ijk - Chr_r_minus_eps$Ijk ) / (2 * EPS);
 Chr$Ijk,a = ( Chr_a_plus_eps$Ijk - Chr_a_minus_eps$Ijk ) / (2 * EPS);
 Chr$Ijk,p = 0.0;
 R$Ijkl = Chr$Ijl,k - Chr$Ijk,l + Chr$Imk * Chr$Mjl - Chr$Iml * Chr$Mjk;
 R$jk = R$Ijik;
 R = g$IL * R$li;
$$
 return( R );
}

main()
{
 int x,y;
 double radius,theta,phi,R;
 double kerr_curvature();
 FILE *pfile;

 pfile = fopen("results.dat","w");
 theta = PI/2;
 array4x4 = matrix(0,3,0,3);
 for (x = -10; x <= 10; x++) {
 for (y = -10; y <= 10; y++) {
    radius = sqrt((double)(x*x) + (double)(y*y));
    phi = atan2((double)(y), (double)(x));
    if (radius > 0) 
      {
       R = kerr_curvature(radius, theta, phi);
       fprintf(pfile,"%d %d %g \n",x,y,R);
      }
    else
      fprintf(pfile,"%d %d %g \n",0,0,0.0);
 } }
 free_matrix(array4x4,0,3,0,3);
 fclose(pfile);
}


