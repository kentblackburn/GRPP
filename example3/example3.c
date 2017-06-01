/* 
   	File:		example3.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Charged particle in Electromagnetic field

	Copyright (c) 1993-1998
	By James Kent Blackburn
	All Rights Reserved
*/


/* Open  GRPP Block */
/* set up coordinate system and tensor indices for GRPP */
/* COORDINATES: TtXxYyZz */
/* DIMENSION: 4 */

#include "example3_grpp.h" 

/* INDICES: IiJjKkLlMmNn */
/* TOTAL: 6 */

/* Close GRPP Block */

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
 
/* Open  GRPP Block */
 TENSOR_C r_C ;
 TENSOR_C u_C ;
 TENSOR_C a_C ;
/* Close GRPP Block */

 dtau = (1.0 / FREQ) / 100000;
 t = x = y = z = 0;
 ux = uy = uz = 0; ut = 1;
 
/* Open  GRPP Block */

 r_C.T = t;
 r_C.X = x;
 r_C.Y = y;
 r_C.Z = z;

 u_C.T = ut;
 u_C.X = ux;
 u_C.Y = uy;
 u_C.Z = uz;

/* Close GRPP Block */
 printf(" %f %f %f %f %f \n",tau, t, x, y, z);
 for (nprint=0; nprint<100; nprint++)
    {
     for (isteps=0; isteps<2000; isteps++)
        {
         acceleration( t, u_C, &a_C );
         
/* Open  GRPP Block */

 u_C.T += ( a_C.T * dtau );
 u_C.X += ( a_C.X * dtau );
 u_C.Y += ( a_C.Y * dtau );
 u_C.Z += ( a_C.Z * dtau );

 r_C.T += ( u_C.T * dtau );
 r_C.X += ( u_C.X * dtau );
 r_C.Y += ( u_C.Y * dtau );
 r_C.Z += ( u_C.Z * dtau );

 t = r_C.T;

/* Close GRPP Block */
         tau += dtau;
        }
     
/* Open  GRPP Block */

 ut = u_C.T;
 ux = u_C.X;
 uy = u_C.Y;
 uz = u_C.Z;

 t = r_C.T;
 x = r_C.X;
 y = r_C.Y;
 z = r_C.Z;

/* Close GRPP Block */
     printf(" %f %f %f %f %f \n",tau, t, x, y, z);
    }
}

void acceleration( t, u_C, a_C )
double t;

/* Open  GRPP Block */
 TENSOR_C u_C ;
 TENSOR_C *a_C ;
/* Close GRPP Block */
{
 double q,m,omega,phi,
        E_o,E_x,E_y,E_z,
        B_o,B_x,B_y,B_z;

/* Open  GRPP Block */

 TENSOR_Cc Faraday_Cc ;

/* Close GRPP Block */

 omega = 2 * PI * FREQ;
 phi = PI / 2;
 q = CHARGE; m = MASS;
 E_o = 10000;
 B_o = 10;
 E_z = 0;
 B_z = 0;

/* Open  GRPP Block */
/* Initialize the tensors */
 Faraday_Cc.Tt = 0;
 Faraday_Cc.Xx = 0;
 Faraday_Cc.Yy = 0;
 Faraday_Cc.Zz = 0;

 Faraday_Cc.Tz = ( Faraday_Cc.Zt = E_z );

 Faraday_Cc.Xy = -( Faraday_Cc.Yx = -B_z );

/* Close GRPP Block */
 E_x =  E_o * cos( omega * t );
 B_y =  B_o * cos( omega * t );
 E_y =  E_o * cos( 2 * omega * t + phi); 
 B_x =  B_o * cos( 2 * omega * t + phi); 

/* Open  GRPP Block */

 Faraday_Cc.Tx = ( Faraday_Cc.Xt = E_x );

 Faraday_Cc.Ty = ( Faraday_Cc.Yt = E_y );

 Faraday_Cc.Xz = -( Faraday_Cc.Zx = B_y );

 Faraday_Cc.Yz = -( Faraday_Cc.Zy = -B_x );

 (*a_C).T = ( ( q / m ) * ( Faraday_Cc.Tt * u_C.T + Faraday_Cc.Tx * u_C.X 
            + Faraday_Cc.Ty * u_C.Y + Faraday_Cc.Tz * u_C.Z ) );
 (*a_C).X = ( ( q / m ) * ( Faraday_Cc.Xt * u_C.T + Faraday_Cc.Xx * u_C.X 
            + Faraday_Cc.Xy * u_C.Y + Faraday_Cc.Xz * u_C.Z ) );
 (*a_C).Y = ( ( q / m ) * ( Faraday_Cc.Yt * u_C.T + Faraday_Cc.Yx * u_C.X 
            + Faraday_Cc.Yy * u_C.Y + Faraday_Cc.Yz * u_C.Z ) );
 (*a_C).Z = ( ( q / m ) * ( Faraday_Cc.Zt * u_C.T + Faraday_Cc.Zx * u_C.X 
            + Faraday_Cc.Zy * u_C.Y + Faraday_Cc.Zz * u_C.Z ) );

/* Close GRPP Block */
}
