/* 
   	File:		example3.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Charged particle in Electromagnetic field

	Copyright (c) 1993-1998
	By James Kent Blackburn
	All Rights Reserved
*/

$$ /* set up coordinate system and tensor indices for GRPP */
coordinates$(t,x,y,z);
indices$(i,j,k,l,m,n);
$$

#include <stdio.h>
#include <math.h>

#define PI	3.14159265359
#define MASS	1.0
#define CHARGE	1.0
#define FREQ	100

main()
{
 int nprint,isteps;
 double t,x,y,z,ut,ux,uy,uz,tau,dtau;
 void acceleration();
 $$ rank1 r$I, u$I, a$I; $$

 dtau = (1.0 / FREQ) / 100000;
 t = x = y = z = 0;
 ux = uy = uz = 0; ut = 1;
 $$
 r$I = {t,x,y,z};
 u$I = {ut,ux,uy,uz};
 $$
 printf(" %f %f %f %f %f \n",tau, t, x, y, z);
 for (nprint=0; nprint<100; nprint++)
    {
     for (isteps=0; isteps<2000; isteps++)
        {
         acceleration( t, u_C, &a_C );
         $$
         u$I += a$I * dtau;
         r$I += u$I * dtau;
         t = r$T;
         $$
         tau += dtau;
        }
     $$
     {ut,ux,uy,uz} = u$I;
     {t,x,y,z} = r$I;
     $$
     printf(" %f %f %f %f %f \n",tau, t, x, y, z);
    }
}

void acceleration( t, u_C, a_C )
double t;
$$ rank1 u$I, *a$I; $$
{
 double q,m,omega,phi,
        E_o,E_x,E_y,E_z,
        B_o,B_x,B_y,B_z;
$$
 rank2 Faraday$Ij;
$$

 omega = 2 * PI * FREQ;
 phi = PI / 2;
 q = CHARGE; m = MASS;
 E_o = 10000;
 B_o = 10;
 E_z = 0;
 B_z = 0;
$$ /* Initialize the tensors */
 Faraday$Ii = 0;
 Faraday$Zt =+=  E_z;
 Faraday$Yx =-= -B_z;
$$
 E_x =  E_o * cos( omega * t );
 B_y =  B_o * cos( omega * t );
 E_y =  E_o * cos( 2 * omega * t + phi); 
 B_x =  B_o * cos( 2 * omega * t + phi); 
$$
 Faraday$Xt =+=  E_x;
 Faraday$Yt =+=  E_y;
 Faraday$Zx =-=  B_y;
 Faraday$Zy =-= -B_x;
 *a$I = (q/m) * Faraday$Ij * u$J;
$$
}
