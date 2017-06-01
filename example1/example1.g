/* 
   	File:		example1.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Two dimensional glacier ice flow

	Copyright (c) 1993 
	By James Kent Blackburn
	All Rights Reserved
*/

$$ /* set up coordinate system and tensor indices for GRPP */
coordinates$(x,y);
indices$(i,j,k);
$$

#include <stdio.h>
#include <math.h>

#define SECONDS_PER_MONTH	 2628000
#define SECONDS_PER_YEAR	31536000
#define TIME_STEP		     100

static char *mon[12] = {"Jan","Feb","Mar","Apr",
                        "May","Jun","Jul","Aug",
                        "Sep","Oct","Nov","Dec"};

main()
{
 int month,seconds;
 double t,x,y,ux,uy,dt = ((double)TIME_STEP / (double)SECONDS_PER_YEAR),
                    x_o = 0.0,
                    y_o = 0.0,
                    ux_o = 0.0,
                    uy_o = 0.05,
                    longitudinal_stretching = 0.0001,
                    side_shearing = 0.1,
                    lateral_extension = 0.0001,
	            flowline_turning = 0.1;

$$ /* Define tensors */
rank1 r$I, u_o$i, u$i, u$I;
rank2 eta$IJ, u$i,j;
$$

$$ /* initialize tensors */
 r$I = {x_o, y_o};
 u_o$i = {ux_o, uy_o};
 u$x,x = longitudinal_stretching;
 u$x,y = side_shearing;
 u$y,x = flowline_turning;
 u$y,y = lateral_extension;
 eta$XY =+= 0.0;
 eta$XX = 1.0;
 eta$YY = 1.0;
$$
 x = x_o; y = y_o; ux = ux_o, uy = uy_o;
 printf("Initial:    position[km](%f,%f) velocity[km/yr](%f,%f) \n",x,y,ux,uy);
 for (month=0; month<12; month++) {
    for (seconds=0; seconds<SECONDS_PER_MONTH; seconds+=TIME_STEP) {
       $$ /* integrate */
       u$i = u_o$i + u$i,j * r$J;
       r$I += eta$IJ * u$j * dt;
       $$
       }
    $$ /* place monthly results in ANSI C variables for printing */
    {x,y} = r$I;
    {ux,uy} = u$i;
    $$
 printf("End of %s: position[km](%f,%f) velocity[km/yr](%f,%f) \n",mon[month],x,y,ux,uy);
 } 
}
