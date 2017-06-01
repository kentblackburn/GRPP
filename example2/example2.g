/* 
   	File:		example2.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Test particle free fall into Schwarzchild black hole

	Copyright (c) 1993-1998 
	By James Kent Blackburn
	All Rights Reserved
*/

$$ /* set up coordinate system and tensor indices for GRPP */
coordinates$(t,r,a,p);
indices$(i,j,k,l,m,n);
$$

#include <stdio.h>
#include <math.h>

main()
{
 int i,j,k,line,nsteps;
 double t,r,a,p,ut,ur,ua,up;
 double sin_a,cos_a,R,M,tau,dtau;
 double g_array[4][4];
 
 $$  /* Define tensors */
 rank1 r$I, u$I, accel$I, diag$i, diag$I;
 rank2 g$ij, g$IJ;
 rank3 g$ij,k, Chris$Ijk;
 $$

 M = 1;
 tau = t = 0;		ut = 1;
 r = 8.0;		ur = 0;
 a = 3.14157/2.0;	ua = 0;
 p = 3.14157/4.0;	up = 0;
 sin_a = sin(a);
 cos_a = cos(a);
 dtau = 1.0e-6;
 printf("Tp: %f Tc: %f radius: %f theta: %f phi: %f\n",tau,t,r,a,p);
$$
 g$ij = 0;
 g$IJ = 0;
 g$ij,k = 0;
$$
 for (line=0; line<26; line++) {
 for (nsteps=0;nsteps<1000000;nsteps++) {
 $$
 r$I = {t,r,a,p};
 u$I = {ut,ur,ua,up};
 diag$i = { -(1.0-2.0*M/r), 1.0/(1.0-2.0*M/r), r*r, r*r*sin_a*sin_a };
 g$ii = diag$i;
 diag$I = {-1.0/(1.0-2.0*M/r),(1.0-2.0*M/r),1.0/(r*r),1.0/(r*r*sin_a*sin_a)};
 g$II = diag$I;
 g$tt,r = - 2.0 * M / ( r * r );
 g$rr,r = 2.0 * M / (( r - 2.0 * M ) * ( r - 2.0 * M ));
 g$aa,r = 2.0 * r;
 g$pp,r = 2.0 * r * sin_a * sin_a;
 g$pp,a = 2.0 * r * r * sin_a * cos_a;
 Chris$Ijk = 0.5 * g$IL * ( g$lj,k + g$lk,j - g$kj,l );
 accel$I = - ( Chris$Ijk * u$J * u$K );
 u$I += accel$I * dtau;
 r$I += u$I * dtau;
 {t,r,a,p} = r$I;
 {ut,ur,ua,up} = u$I;
$$
 tau += dtau;
 if (r<=(2*M))
   {
    printf("Passed through event horizon at proper time: %f \n",tau);
    printf("                     and at coordinate time: %f \n",t);
    exit(0);
   }
 }
 printf("Tp: %f Tc: %f radius: %f theta: %f phi: %f\n",tau,t,r,a,p); }
}
