/* 
   	File:		example2.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Test particle free fall into Schwarzchild black hole

	Copyright (c) 1993-1998 
	By James Kent Blackburn
	All Rights Reserved
*/


/* Open  GRPP Block */
/* set up coordinate system and tensor indices for GRPP */
/* COORDINATES: TtRrAaPp */
/* DIMENSION: 4 */

#include "example2_grpp.h" 

/* INDICES: IiJjKkLlMmNn */
/* TOTAL: 6 */

/* Close GRPP Block */

#include <stdio.h>
#include <math.h>

main()
{
 int i,j,k,line,nsteps;
 double t,r,a,p,ut,ur,ua,up;
 double sin_a,cos_a,R,M,tau,dtau;
 double g_array[4][4];
 
 
/* Open  GRPP Block */
/* Define tensors */
 TENSOR_C r_C ;
 TENSOR_C u_C ;
 TENSOR_C accel_C ;
 TENSOR_c diag_c ;
 TENSOR_C diag_C ;

 TENSOR_cc g_cc ;
 TENSOR_CC g_CC ;

 TENSOR_ccc _d_g_ccc ;
 TENSOR_Ccc Chris_Ccc ;

/* Close GRPP Block */

 M = 1;
 tau = t = 0;		ut = 1;
 r = 8.0;		ur = 0;
 a = 3.14157/2.0;	ua = 0;
 p = 3.14157/4.0;	up = 0;
 sin_a = sin(a);
 cos_a = cos(a);
 dtau = 1.0e-6;
 printf("Tp: %f Tc: %f radius: %f theta: %f phi: %f\n",tau,t,r,a,p);

/* Open  GRPP Block */

 g_cc.tt = 0;
 g_cc.tr = 0;
 g_cc.ta = 0;
 g_cc.tp = 0;
 g_cc.rt = 0;
 g_cc.rr = 0;
 g_cc.ra = 0;
 g_cc.rp = 0;
 g_cc.at = 0;
 g_cc.ar = 0;
 g_cc.aa = 0;
 g_cc.ap = 0;
 g_cc.pt = 0;
 g_cc.pr = 0;
 g_cc.pa = 0;
 g_cc.pp = 0;

 g_CC.TT = 0;
 g_CC.TR = 0;
 g_CC.TA = 0;
 g_CC.TP = 0;
 g_CC.RT = 0;
 g_CC.RR = 0;
 g_CC.RA = 0;
 g_CC.RP = 0;
 g_CC.AT = 0;
 g_CC.AR = 0;
 g_CC.AA = 0;
 g_CC.AP = 0;
 g_CC.PT = 0;
 g_CC.PR = 0;
 g_CC.PA = 0;
 g_CC.PP = 0;

 _d_g_ccc.ttt = 0;
 _d_g_ccc.ttr = 0;
 _d_g_ccc.tta = 0;
 _d_g_ccc.ttp = 0;
 _d_g_ccc.trt = 0;
 _d_g_ccc.trr = 0;
 _d_g_ccc.tra = 0;
 _d_g_ccc.trp = 0;
 _d_g_ccc.tat = 0;
 _d_g_ccc.tar = 0;
 _d_g_ccc.taa = 0;
 _d_g_ccc.tap = 0;
 _d_g_ccc.tpt = 0;
 _d_g_ccc.tpr = 0;
 _d_g_ccc.tpa = 0;
 _d_g_ccc.tpp = 0;
 _d_g_ccc.rtt = 0;
 _d_g_ccc.rtr = 0;
 _d_g_ccc.rta = 0;
 _d_g_ccc.rtp = 0;
 _d_g_ccc.rrt = 0;
 _d_g_ccc.rrr = 0;
 _d_g_ccc.rra = 0;
 _d_g_ccc.rrp = 0;
 _d_g_ccc.rat = 0;
 _d_g_ccc.rar = 0;
 _d_g_ccc.raa = 0;
 _d_g_ccc.rap = 0;
 _d_g_ccc.rpt = 0;
 _d_g_ccc.rpr = 0;
 _d_g_ccc.rpa = 0;
 _d_g_ccc.rpp = 0;
 _d_g_ccc.att = 0;
 _d_g_ccc.atr = 0;
 _d_g_ccc.ata = 0;
 _d_g_ccc.atp = 0;
 _d_g_ccc.art = 0;
 _d_g_ccc.arr = 0;
 _d_g_ccc.ara = 0;
 _d_g_ccc.arp = 0;
 _d_g_ccc.aat = 0;
 _d_g_ccc.aar = 0;
 _d_g_ccc.aaa = 0;
 _d_g_ccc.aap = 0;
 _d_g_ccc.apt = 0;
 _d_g_ccc.apr = 0;
 _d_g_ccc.apa = 0;
 _d_g_ccc.app = 0;
 _d_g_ccc.ptt = 0;
 _d_g_ccc.ptr = 0;
 _d_g_ccc.pta = 0;
 _d_g_ccc.ptp = 0;
 _d_g_ccc.prt = 0;
 _d_g_ccc.prr = 0;
 _d_g_ccc.pra = 0;
 _d_g_ccc.prp = 0;
 _d_g_ccc.pat = 0;
 _d_g_ccc.par = 0;
 _d_g_ccc.paa = 0;
 _d_g_ccc.pap = 0;
 _d_g_ccc.ppt = 0;
 _d_g_ccc.ppr = 0;
 _d_g_ccc.ppa = 0;
 _d_g_ccc.ppp = 0;

/* Close GRPP Block */
 for (line=0; line<26; line++) {
 for (nsteps=0;nsteps<1000000;nsteps++) {
 
/* Open  GRPP Block */

 r_C.T = t;
 r_C.R = r;
 r_C.A = a;
 r_C.P = p;

 u_C.T = ut;
 u_C.R = ur;
 u_C.A = ua;
 u_C.P = up;

 diag_c.t =  -(1.0-2.0*M/r);
 diag_c.r =  1.0/(1.0-2.0*M/r);
 diag_c.a =  r*r;
 diag_c.p =  r*r*sin_a*sin_a ;

 g_cc.tt = diag_c.t;
 g_cc.rr = diag_c.r;
 g_cc.aa = diag_c.a;
 g_cc.pp = diag_c.p;

 diag_C.T = -1.0/(1.0-2.0*M/r);
 diag_C.R = (1.0-2.0*M/r);
 diag_C.A = 1.0/(r*r);
 diag_C.P = 1.0/(r*r*sin_a*sin_a);

 g_CC.TT = diag_C.T;
 g_CC.RR = diag_C.R;
 g_CC.AA = diag_C.A;
 g_CC.PP = diag_C.P;

 _d_g_ccc.ttr = ( -2.0 * M / ( r * r ) );

 _d_g_ccc.rrr = ( 2.0 * M / ( ( r - ( 2.0 * M ) ) * ( r - ( 2.0 * M ) ) ) );

 _d_g_ccc.aar = ( 2.0 * r );

 _d_g_ccc.ppr = ( 2.0 * ( r * ( sin_a * sin_a ) ) );

 _d_g_ccc.ppa = ( 2.0 * ( r * ( r * ( sin_a * cos_a ) ) ) );

 Chris_Ccc.Ttt = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                 - _d_g_ccc.ttt ) + g_CC.TR * ( _d_g_ccc.rtt + _d_g_ccc.rtt 
                 - _d_g_ccc.ttr ) + g_CC.TA * ( _d_g_ccc.att + _d_g_ccc.att 
                 - _d_g_ccc.tta ) + g_CC.TP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                 - _d_g_ccc.ttp ) ) );
 Chris_Ccc.Ttr = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                 - _d_g_ccc.rtt ) + g_CC.TR * ( _d_g_ccc.rtr + _d_g_ccc.rrt 
                 - _d_g_ccc.rtr ) + g_CC.TA * ( _d_g_ccc.atr + _d_g_ccc.art 
                 - _d_g_ccc.rta ) + g_CC.TP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                 - _d_g_ccc.rtp ) ) );
 Chris_Ccc.Tta = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                 - _d_g_ccc.att ) + g_CC.TR * ( _d_g_ccc.rta + _d_g_ccc.rat 
                 - _d_g_ccc.atr ) + g_CC.TA * ( _d_g_ccc.ata + _d_g_ccc.aat 
                 - _d_g_ccc.ata ) + g_CC.TP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                 - _d_g_ccc.atp ) ) );
 Chris_Ccc.Ttp = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                 - _d_g_ccc.ptt ) + g_CC.TR * ( _d_g_ccc.rtp + _d_g_ccc.rpt 
                 - _d_g_ccc.ptr ) + g_CC.TA * ( _d_g_ccc.atp + _d_g_ccc.apt 
                 - _d_g_ccc.pta ) + g_CC.TP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                 - _d_g_ccc.ptp ) ) );
 Chris_Ccc.Trt = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                 - _d_g_ccc.trt ) + g_CC.TR * ( _d_g_ccc.rrt + _d_g_ccc.rtr 
                 - _d_g_ccc.trr ) + g_CC.TA * ( _d_g_ccc.art + _d_g_ccc.atr 
                 - _d_g_ccc.tra ) + g_CC.TP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                 - _d_g_ccc.trp ) ) );
 Chris_Ccc.Trr = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                 - _d_g_ccc.rrt ) + g_CC.TR * ( _d_g_ccc.rrr + _d_g_ccc.rrr 
                 - _d_g_ccc.rrr ) + g_CC.TA * ( _d_g_ccc.arr + _d_g_ccc.arr 
                 - _d_g_ccc.rra ) + g_CC.TP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                 - _d_g_ccc.rrp ) ) );
 Chris_Ccc.Tra = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                 - _d_g_ccc.art ) + g_CC.TR * ( _d_g_ccc.rra + _d_g_ccc.rar 
                 - _d_g_ccc.arr ) + g_CC.TA * ( _d_g_ccc.ara + _d_g_ccc.aar 
                 - _d_g_ccc.ara ) + g_CC.TP * ( _d_g_ccc.pra + _d_g_ccc.par 
                 - _d_g_ccc.arp ) ) );
 Chris_Ccc.Trp = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                 - _d_g_ccc.prt ) + g_CC.TR * ( _d_g_ccc.rrp + _d_g_ccc.rpr 
                 - _d_g_ccc.prr ) + g_CC.TA * ( _d_g_ccc.arp + _d_g_ccc.apr 
                 - _d_g_ccc.pra ) + g_CC.TP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                 - _d_g_ccc.prp ) ) );
 Chris_Ccc.Tat = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                 - _d_g_ccc.tat ) + g_CC.TR * ( _d_g_ccc.rat + _d_g_ccc.rta 
                 - _d_g_ccc.tar ) + g_CC.TA * ( _d_g_ccc.aat + _d_g_ccc.ata 
                 - _d_g_ccc.taa ) + g_CC.TP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                 - _d_g_ccc.tap ) ) );
 Chris_Ccc.Tar = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                 - _d_g_ccc.rat ) + g_CC.TR * ( _d_g_ccc.rar + _d_g_ccc.rra 
                 - _d_g_ccc.rar ) + g_CC.TA * ( _d_g_ccc.aar + _d_g_ccc.ara 
                 - _d_g_ccc.raa ) + g_CC.TP * ( _d_g_ccc.par + _d_g_ccc.pra 
                 - _d_g_ccc.rap ) ) );
 Chris_Ccc.Taa = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                 - _d_g_ccc.aat ) + g_CC.TR * ( _d_g_ccc.raa + _d_g_ccc.raa 
                 - _d_g_ccc.aar ) + g_CC.TA * ( _d_g_ccc.aaa + _d_g_ccc.aaa 
                 - _d_g_ccc.aaa ) + g_CC.TP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                 - _d_g_ccc.aap ) ) );
 Chris_Ccc.Tap = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                 - _d_g_ccc.pat ) + g_CC.TR * ( _d_g_ccc.rap + _d_g_ccc.rpa 
                 - _d_g_ccc.par ) + g_CC.TA * ( _d_g_ccc.aap + _d_g_ccc.apa 
                 - _d_g_ccc.paa ) + g_CC.TP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                 - _d_g_ccc.pap ) ) );
 Chris_Ccc.Tpt = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                 - _d_g_ccc.tpt ) + g_CC.TR * ( _d_g_ccc.rpt + _d_g_ccc.rtp 
                 - _d_g_ccc.tpr ) + g_CC.TA * ( _d_g_ccc.apt + _d_g_ccc.atp 
                 - _d_g_ccc.tpa ) + g_CC.TP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                 - _d_g_ccc.tpp ) ) );
 Chris_Ccc.Tpr = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                 - _d_g_ccc.rpt ) + g_CC.TR * ( _d_g_ccc.rpr + _d_g_ccc.rrp 
                 - _d_g_ccc.rpr ) + g_CC.TA * ( _d_g_ccc.apr + _d_g_ccc.arp 
                 - _d_g_ccc.rpa ) + g_CC.TP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                 - _d_g_ccc.rpp ) ) );
 Chris_Ccc.Tpa = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                 - _d_g_ccc.apt ) + g_CC.TR * ( _d_g_ccc.rpa + _d_g_ccc.rap 
                 - _d_g_ccc.apr ) + g_CC.TA * ( _d_g_ccc.apa + _d_g_ccc.aap 
                 - _d_g_ccc.apa ) + g_CC.TP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                 - _d_g_ccc.app ) ) );
 Chris_Ccc.Tpp = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                 - _d_g_ccc.ppt ) + g_CC.TR * ( _d_g_ccc.rpp + _d_g_ccc.rpp 
                 - _d_g_ccc.ppr ) + g_CC.TA * ( _d_g_ccc.app + _d_g_ccc.app 
                 - _d_g_ccc.ppa ) + g_CC.TP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                 - _d_g_ccc.ppp ) ) );
 Chris_Ccc.Rtt = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                 - _d_g_ccc.ttt ) + g_CC.RR * ( _d_g_ccc.rtt + _d_g_ccc.rtt 
                 - _d_g_ccc.ttr ) + g_CC.RA * ( _d_g_ccc.att + _d_g_ccc.att 
                 - _d_g_ccc.tta ) + g_CC.RP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                 - _d_g_ccc.ttp ) ) );
 Chris_Ccc.Rtr = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                 - _d_g_ccc.rtt ) + g_CC.RR * ( _d_g_ccc.rtr + _d_g_ccc.rrt 
                 - _d_g_ccc.rtr ) + g_CC.RA * ( _d_g_ccc.atr + _d_g_ccc.art 
                 - _d_g_ccc.rta ) + g_CC.RP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                 - _d_g_ccc.rtp ) ) );
 Chris_Ccc.Rta = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                 - _d_g_ccc.att ) + g_CC.RR * ( _d_g_ccc.rta + _d_g_ccc.rat 
                 - _d_g_ccc.atr ) + g_CC.RA * ( _d_g_ccc.ata + _d_g_ccc.aat 
                 - _d_g_ccc.ata ) + g_CC.RP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                 - _d_g_ccc.atp ) ) );
 Chris_Ccc.Rtp = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                 - _d_g_ccc.ptt ) + g_CC.RR * ( _d_g_ccc.rtp + _d_g_ccc.rpt 
                 - _d_g_ccc.ptr ) + g_CC.RA * ( _d_g_ccc.atp + _d_g_ccc.apt 
                 - _d_g_ccc.pta ) + g_CC.RP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                 - _d_g_ccc.ptp ) ) );
 Chris_Ccc.Rrt = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                 - _d_g_ccc.trt ) + g_CC.RR * ( _d_g_ccc.rrt + _d_g_ccc.rtr 
                 - _d_g_ccc.trr ) + g_CC.RA * ( _d_g_ccc.art + _d_g_ccc.atr 
                 - _d_g_ccc.tra ) + g_CC.RP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                 - _d_g_ccc.trp ) ) );
 Chris_Ccc.Rrr = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                 - _d_g_ccc.rrt ) + g_CC.RR * ( _d_g_ccc.rrr + _d_g_ccc.rrr 
                 - _d_g_ccc.rrr ) + g_CC.RA * ( _d_g_ccc.arr + _d_g_ccc.arr 
                 - _d_g_ccc.rra ) + g_CC.RP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                 - _d_g_ccc.rrp ) ) );
 Chris_Ccc.Rra = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                 - _d_g_ccc.art ) + g_CC.RR * ( _d_g_ccc.rra + _d_g_ccc.rar 
                 - _d_g_ccc.arr ) + g_CC.RA * ( _d_g_ccc.ara + _d_g_ccc.aar 
                 - _d_g_ccc.ara ) + g_CC.RP * ( _d_g_ccc.pra + _d_g_ccc.par 
                 - _d_g_ccc.arp ) ) );
 Chris_Ccc.Rrp = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                 - _d_g_ccc.prt ) + g_CC.RR * ( _d_g_ccc.rrp + _d_g_ccc.rpr 
                 - _d_g_ccc.prr ) + g_CC.RA * ( _d_g_ccc.arp + _d_g_ccc.apr 
                 - _d_g_ccc.pra ) + g_CC.RP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                 - _d_g_ccc.prp ) ) );
 Chris_Ccc.Rat = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                 - _d_g_ccc.tat ) + g_CC.RR * ( _d_g_ccc.rat + _d_g_ccc.rta 
                 - _d_g_ccc.tar ) + g_CC.RA * ( _d_g_ccc.aat + _d_g_ccc.ata 
                 - _d_g_ccc.taa ) + g_CC.RP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                 - _d_g_ccc.tap ) ) );
 Chris_Ccc.Rar = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                 - _d_g_ccc.rat ) + g_CC.RR * ( _d_g_ccc.rar + _d_g_ccc.rra 
                 - _d_g_ccc.rar ) + g_CC.RA * ( _d_g_ccc.aar + _d_g_ccc.ara 
                 - _d_g_ccc.raa ) + g_CC.RP * ( _d_g_ccc.par + _d_g_ccc.pra 
                 - _d_g_ccc.rap ) ) );
 Chris_Ccc.Raa = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                 - _d_g_ccc.aat ) + g_CC.RR * ( _d_g_ccc.raa + _d_g_ccc.raa 
                 - _d_g_ccc.aar ) + g_CC.RA * ( _d_g_ccc.aaa + _d_g_ccc.aaa 
                 - _d_g_ccc.aaa ) + g_CC.RP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                 - _d_g_ccc.aap ) ) );
 Chris_Ccc.Rap = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                 - _d_g_ccc.pat ) + g_CC.RR * ( _d_g_ccc.rap + _d_g_ccc.rpa 
                 - _d_g_ccc.par ) + g_CC.RA * ( _d_g_ccc.aap + _d_g_ccc.apa 
                 - _d_g_ccc.paa ) + g_CC.RP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                 - _d_g_ccc.pap ) ) );
 Chris_Ccc.Rpt = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                 - _d_g_ccc.tpt ) + g_CC.RR * ( _d_g_ccc.rpt + _d_g_ccc.rtp 
                 - _d_g_ccc.tpr ) + g_CC.RA * ( _d_g_ccc.apt + _d_g_ccc.atp 
                 - _d_g_ccc.tpa ) + g_CC.RP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                 - _d_g_ccc.tpp ) ) );
 Chris_Ccc.Rpr = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                 - _d_g_ccc.rpt ) + g_CC.RR * ( _d_g_ccc.rpr + _d_g_ccc.rrp 
                 - _d_g_ccc.rpr ) + g_CC.RA * ( _d_g_ccc.apr + _d_g_ccc.arp 
                 - _d_g_ccc.rpa ) + g_CC.RP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                 - _d_g_ccc.rpp ) ) );
 Chris_Ccc.Rpa = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                 - _d_g_ccc.apt ) + g_CC.RR * ( _d_g_ccc.rpa + _d_g_ccc.rap 
                 - _d_g_ccc.apr ) + g_CC.RA * ( _d_g_ccc.apa + _d_g_ccc.aap 
                 - _d_g_ccc.apa ) + g_CC.RP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                 - _d_g_ccc.app ) ) );
 Chris_Ccc.Rpp = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                 - _d_g_ccc.ppt ) + g_CC.RR * ( _d_g_ccc.rpp + _d_g_ccc.rpp 
                 - _d_g_ccc.ppr ) + g_CC.RA * ( _d_g_ccc.app + _d_g_ccc.app 
                 - _d_g_ccc.ppa ) + g_CC.RP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                 - _d_g_ccc.ppp ) ) );
 Chris_Ccc.Att = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                 - _d_g_ccc.ttt ) + g_CC.AR * ( _d_g_ccc.rtt + _d_g_ccc.rtt 
                 - _d_g_ccc.ttr ) + g_CC.AA * ( _d_g_ccc.att + _d_g_ccc.att 
                 - _d_g_ccc.tta ) + g_CC.AP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                 - _d_g_ccc.ttp ) ) );
 Chris_Ccc.Atr = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                 - _d_g_ccc.rtt ) + g_CC.AR * ( _d_g_ccc.rtr + _d_g_ccc.rrt 
                 - _d_g_ccc.rtr ) + g_CC.AA * ( _d_g_ccc.atr + _d_g_ccc.art 
                 - _d_g_ccc.rta ) + g_CC.AP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                 - _d_g_ccc.rtp ) ) );
 Chris_Ccc.Ata = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                 - _d_g_ccc.att ) + g_CC.AR * ( _d_g_ccc.rta + _d_g_ccc.rat 
                 - _d_g_ccc.atr ) + g_CC.AA * ( _d_g_ccc.ata + _d_g_ccc.aat 
                 - _d_g_ccc.ata ) + g_CC.AP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                 - _d_g_ccc.atp ) ) );
 Chris_Ccc.Atp = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                 - _d_g_ccc.ptt ) + g_CC.AR * ( _d_g_ccc.rtp + _d_g_ccc.rpt 
                 - _d_g_ccc.ptr ) + g_CC.AA * ( _d_g_ccc.atp + _d_g_ccc.apt 
                 - _d_g_ccc.pta ) + g_CC.AP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                 - _d_g_ccc.ptp ) ) );
 Chris_Ccc.Art = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                 - _d_g_ccc.trt ) + g_CC.AR * ( _d_g_ccc.rrt + _d_g_ccc.rtr 
                 - _d_g_ccc.trr ) + g_CC.AA * ( _d_g_ccc.art + _d_g_ccc.atr 
                 - _d_g_ccc.tra ) + g_CC.AP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                 - _d_g_ccc.trp ) ) );
 Chris_Ccc.Arr = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                 - _d_g_ccc.rrt ) + g_CC.AR * ( _d_g_ccc.rrr + _d_g_ccc.rrr 
                 - _d_g_ccc.rrr ) + g_CC.AA * ( _d_g_ccc.arr + _d_g_ccc.arr 
                 - _d_g_ccc.rra ) + g_CC.AP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                 - _d_g_ccc.rrp ) ) );
 Chris_Ccc.Ara = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                 - _d_g_ccc.art ) + g_CC.AR * ( _d_g_ccc.rra + _d_g_ccc.rar 
                 - _d_g_ccc.arr ) + g_CC.AA * ( _d_g_ccc.ara + _d_g_ccc.aar 
                 - _d_g_ccc.ara ) + g_CC.AP * ( _d_g_ccc.pra + _d_g_ccc.par 
                 - _d_g_ccc.arp ) ) );
 Chris_Ccc.Arp = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                 - _d_g_ccc.prt ) + g_CC.AR * ( _d_g_ccc.rrp + _d_g_ccc.rpr 
                 - _d_g_ccc.prr ) + g_CC.AA * ( _d_g_ccc.arp + _d_g_ccc.apr 
                 - _d_g_ccc.pra ) + g_CC.AP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                 - _d_g_ccc.prp ) ) );
 Chris_Ccc.Aat = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                 - _d_g_ccc.tat ) + g_CC.AR * ( _d_g_ccc.rat + _d_g_ccc.rta 
                 - _d_g_ccc.tar ) + g_CC.AA * ( _d_g_ccc.aat + _d_g_ccc.ata 
                 - _d_g_ccc.taa ) + g_CC.AP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                 - _d_g_ccc.tap ) ) );
 Chris_Ccc.Aar = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                 - _d_g_ccc.rat ) + g_CC.AR * ( _d_g_ccc.rar + _d_g_ccc.rra 
                 - _d_g_ccc.rar ) + g_CC.AA * ( _d_g_ccc.aar + _d_g_ccc.ara 
                 - _d_g_ccc.raa ) + g_CC.AP * ( _d_g_ccc.par + _d_g_ccc.pra 
                 - _d_g_ccc.rap ) ) );
 Chris_Ccc.Aaa = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                 - _d_g_ccc.aat ) + g_CC.AR * ( _d_g_ccc.raa + _d_g_ccc.raa 
                 - _d_g_ccc.aar ) + g_CC.AA * ( _d_g_ccc.aaa + _d_g_ccc.aaa 
                 - _d_g_ccc.aaa ) + g_CC.AP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                 - _d_g_ccc.aap ) ) );
 Chris_Ccc.Aap = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                 - _d_g_ccc.pat ) + g_CC.AR * ( _d_g_ccc.rap + _d_g_ccc.rpa 
                 - _d_g_ccc.par ) + g_CC.AA * ( _d_g_ccc.aap + _d_g_ccc.apa 
                 - _d_g_ccc.paa ) + g_CC.AP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                 - _d_g_ccc.pap ) ) );
 Chris_Ccc.Apt = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                 - _d_g_ccc.tpt ) + g_CC.AR * ( _d_g_ccc.rpt + _d_g_ccc.rtp 
                 - _d_g_ccc.tpr ) + g_CC.AA * ( _d_g_ccc.apt + _d_g_ccc.atp 
                 - _d_g_ccc.tpa ) + g_CC.AP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                 - _d_g_ccc.tpp ) ) );
 Chris_Ccc.Apr = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                 - _d_g_ccc.rpt ) + g_CC.AR * ( _d_g_ccc.rpr + _d_g_ccc.rrp 
                 - _d_g_ccc.rpr ) + g_CC.AA * ( _d_g_ccc.apr + _d_g_ccc.arp 
                 - _d_g_ccc.rpa ) + g_CC.AP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                 - _d_g_ccc.rpp ) ) );
 Chris_Ccc.Apa = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                 - _d_g_ccc.apt ) + g_CC.AR * ( _d_g_ccc.rpa + _d_g_ccc.rap 
                 - _d_g_ccc.apr ) + g_CC.AA * ( _d_g_ccc.apa + _d_g_ccc.aap 
                 - _d_g_ccc.apa ) + g_CC.AP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                 - _d_g_ccc.app ) ) );
 Chris_Ccc.App = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                 - _d_g_ccc.ppt ) + g_CC.AR * ( _d_g_ccc.rpp + _d_g_ccc.rpp 
                 - _d_g_ccc.ppr ) + g_CC.AA * ( _d_g_ccc.app + _d_g_ccc.app 
                 - _d_g_ccc.ppa ) + g_CC.AP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                 - _d_g_ccc.ppp ) ) );
 Chris_Ccc.Ptt = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                 - _d_g_ccc.ttt ) + g_CC.PR * ( _d_g_ccc.rtt + _d_g_ccc.rtt 
                 - _d_g_ccc.ttr ) + g_CC.PA * ( _d_g_ccc.att + _d_g_ccc.att 
                 - _d_g_ccc.tta ) + g_CC.PP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                 - _d_g_ccc.ttp ) ) );
 Chris_Ccc.Ptr = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                 - _d_g_ccc.rtt ) + g_CC.PR * ( _d_g_ccc.rtr + _d_g_ccc.rrt 
                 - _d_g_ccc.rtr ) + g_CC.PA * ( _d_g_ccc.atr + _d_g_ccc.art 
                 - _d_g_ccc.rta ) + g_CC.PP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                 - _d_g_ccc.rtp ) ) );
 Chris_Ccc.Pta = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                 - _d_g_ccc.att ) + g_CC.PR * ( _d_g_ccc.rta + _d_g_ccc.rat 
                 - _d_g_ccc.atr ) + g_CC.PA * ( _d_g_ccc.ata + _d_g_ccc.aat 
                 - _d_g_ccc.ata ) + g_CC.PP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                 - _d_g_ccc.atp ) ) );
 Chris_Ccc.Ptp = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                 - _d_g_ccc.ptt ) + g_CC.PR * ( _d_g_ccc.rtp + _d_g_ccc.rpt 
                 - _d_g_ccc.ptr ) + g_CC.PA * ( _d_g_ccc.atp + _d_g_ccc.apt 
                 - _d_g_ccc.pta ) + g_CC.PP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                 - _d_g_ccc.ptp ) ) );
 Chris_Ccc.Prt = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                 - _d_g_ccc.trt ) + g_CC.PR * ( _d_g_ccc.rrt + _d_g_ccc.rtr 
                 - _d_g_ccc.trr ) + g_CC.PA * ( _d_g_ccc.art + _d_g_ccc.atr 
                 - _d_g_ccc.tra ) + g_CC.PP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                 - _d_g_ccc.trp ) ) );
 Chris_Ccc.Prr = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                 - _d_g_ccc.rrt ) + g_CC.PR * ( _d_g_ccc.rrr + _d_g_ccc.rrr 
                 - _d_g_ccc.rrr ) + g_CC.PA * ( _d_g_ccc.arr + _d_g_ccc.arr 
                 - _d_g_ccc.rra ) + g_CC.PP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                 - _d_g_ccc.rrp ) ) );
 Chris_Ccc.Pra = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                 - _d_g_ccc.art ) + g_CC.PR * ( _d_g_ccc.rra + _d_g_ccc.rar 
                 - _d_g_ccc.arr ) + g_CC.PA * ( _d_g_ccc.ara + _d_g_ccc.aar 
                 - _d_g_ccc.ara ) + g_CC.PP * ( _d_g_ccc.pra + _d_g_ccc.par 
                 - _d_g_ccc.arp ) ) );
 Chris_Ccc.Prp = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                 - _d_g_ccc.prt ) + g_CC.PR * ( _d_g_ccc.rrp + _d_g_ccc.rpr 
                 - _d_g_ccc.prr ) + g_CC.PA * ( _d_g_ccc.arp + _d_g_ccc.apr 
                 - _d_g_ccc.pra ) + g_CC.PP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                 - _d_g_ccc.prp ) ) );
 Chris_Ccc.Pat = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                 - _d_g_ccc.tat ) + g_CC.PR * ( _d_g_ccc.rat + _d_g_ccc.rta 
                 - _d_g_ccc.tar ) + g_CC.PA * ( _d_g_ccc.aat + _d_g_ccc.ata 
                 - _d_g_ccc.taa ) + g_CC.PP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                 - _d_g_ccc.tap ) ) );
 Chris_Ccc.Par = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                 - _d_g_ccc.rat ) + g_CC.PR * ( _d_g_ccc.rar + _d_g_ccc.rra 
                 - _d_g_ccc.rar ) + g_CC.PA * ( _d_g_ccc.aar + _d_g_ccc.ara 
                 - _d_g_ccc.raa ) + g_CC.PP * ( _d_g_ccc.par + _d_g_ccc.pra 
                 - _d_g_ccc.rap ) ) );
 Chris_Ccc.Paa = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                 - _d_g_ccc.aat ) + g_CC.PR * ( _d_g_ccc.raa + _d_g_ccc.raa 
                 - _d_g_ccc.aar ) + g_CC.PA * ( _d_g_ccc.aaa + _d_g_ccc.aaa 
                 - _d_g_ccc.aaa ) + g_CC.PP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                 - _d_g_ccc.aap ) ) );
 Chris_Ccc.Pap = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                 - _d_g_ccc.pat ) + g_CC.PR * ( _d_g_ccc.rap + _d_g_ccc.rpa 
                 - _d_g_ccc.par ) + g_CC.PA * ( _d_g_ccc.aap + _d_g_ccc.apa 
                 - _d_g_ccc.paa ) + g_CC.PP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                 - _d_g_ccc.pap ) ) );
 Chris_Ccc.Ppt = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                 - _d_g_ccc.tpt ) + g_CC.PR * ( _d_g_ccc.rpt + _d_g_ccc.rtp 
                 - _d_g_ccc.tpr ) + g_CC.PA * ( _d_g_ccc.apt + _d_g_ccc.atp 
                 - _d_g_ccc.tpa ) + g_CC.PP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                 - _d_g_ccc.tpp ) ) );
 Chris_Ccc.Ppr = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                 - _d_g_ccc.rpt ) + g_CC.PR * ( _d_g_ccc.rpr + _d_g_ccc.rrp 
                 - _d_g_ccc.rpr ) + g_CC.PA * ( _d_g_ccc.apr + _d_g_ccc.arp 
                 - _d_g_ccc.rpa ) + g_CC.PP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                 - _d_g_ccc.rpp ) ) );
 Chris_Ccc.Ppa = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                 - _d_g_ccc.apt ) + g_CC.PR * ( _d_g_ccc.rpa + _d_g_ccc.rap 
                 - _d_g_ccc.apr ) + g_CC.PA * ( _d_g_ccc.apa + _d_g_ccc.aap 
                 - _d_g_ccc.apa ) + g_CC.PP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                 - _d_g_ccc.app ) ) );
 Chris_Ccc.Ppp = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                 - _d_g_ccc.ppt ) + g_CC.PR * ( _d_g_ccc.rpp + _d_g_ccc.rpp 
                 - _d_g_ccc.ppr ) + g_CC.PA * ( _d_g_ccc.app + _d_g_ccc.app 
                 - _d_g_ccc.ppa ) + g_CC.PP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                 - _d_g_ccc.ppp ) ) );

 accel_C.T = -( Chris_Ccc.Ttt * ( u_C.T * u_C.T ) + Chris_Ccc.Ttr * ( u_C.T 
             * u_C.R ) + Chris_Ccc.Tta * ( u_C.T * u_C.A ) + Chris_Ccc.Ttp 
             * ( u_C.T * u_C.P ) + Chris_Ccc.Trt * ( u_C.R * u_C.T ) 
             + Chris_Ccc.Trr * ( u_C.R * u_C.R ) + Chris_Ccc.Tra * ( u_C.R 
             * u_C.A ) + Chris_Ccc.Trp * ( u_C.R * u_C.P ) + Chris_Ccc.Tat 
             * ( u_C.A * u_C.T ) + Chris_Ccc.Tar * ( u_C.A * u_C.R ) 
             + Chris_Ccc.Taa * ( u_C.A * u_C.A ) + Chris_Ccc.Tap * ( u_C.A 
             * u_C.P ) + Chris_Ccc.Tpt * ( u_C.P * u_C.T ) + Chris_Ccc.Tpr 
             * ( u_C.P * u_C.R ) + Chris_Ccc.Tpa * ( u_C.P * u_C.A ) 
             + Chris_Ccc.Tpp * ( u_C.P * u_C.P ) );
 accel_C.R = -( Chris_Ccc.Rtt * ( u_C.T * u_C.T ) + Chris_Ccc.Rtr * ( u_C.T 
             * u_C.R ) + Chris_Ccc.Rta * ( u_C.T * u_C.A ) + Chris_Ccc.Rtp 
             * ( u_C.T * u_C.P ) + Chris_Ccc.Rrt * ( u_C.R * u_C.T ) 
             + Chris_Ccc.Rrr * ( u_C.R * u_C.R ) + Chris_Ccc.Rra * ( u_C.R 
             * u_C.A ) + Chris_Ccc.Rrp * ( u_C.R * u_C.P ) + Chris_Ccc.Rat 
             * ( u_C.A * u_C.T ) + Chris_Ccc.Rar * ( u_C.A * u_C.R ) 
             + Chris_Ccc.Raa * ( u_C.A * u_C.A ) + Chris_Ccc.Rap * ( u_C.A 
             * u_C.P ) + Chris_Ccc.Rpt * ( u_C.P * u_C.T ) + Chris_Ccc.Rpr 
             * ( u_C.P * u_C.R ) + Chris_Ccc.Rpa * ( u_C.P * u_C.A ) 
             + Chris_Ccc.Rpp * ( u_C.P * u_C.P ) );
 accel_C.A = -( Chris_Ccc.Att * ( u_C.T * u_C.T ) + Chris_Ccc.Atr * ( u_C.T 
             * u_C.R ) + Chris_Ccc.Ata * ( u_C.T * u_C.A ) + Chris_Ccc.Atp 
             * ( u_C.T * u_C.P ) + Chris_Ccc.Art * ( u_C.R * u_C.T ) 
             + Chris_Ccc.Arr * ( u_C.R * u_C.R ) + Chris_Ccc.Ara * ( u_C.R 
             * u_C.A ) + Chris_Ccc.Arp * ( u_C.R * u_C.P ) + Chris_Ccc.Aat 
             * ( u_C.A * u_C.T ) + Chris_Ccc.Aar * ( u_C.A * u_C.R ) 
             + Chris_Ccc.Aaa * ( u_C.A * u_C.A ) + Chris_Ccc.Aap * ( u_C.A 
             * u_C.P ) + Chris_Ccc.Apt * ( u_C.P * u_C.T ) + Chris_Ccc.Apr 
             * ( u_C.P * u_C.R ) + Chris_Ccc.Apa * ( u_C.P * u_C.A ) 
             + Chris_Ccc.App * ( u_C.P * u_C.P ) );
 accel_C.P = -( Chris_Ccc.Ptt * ( u_C.T * u_C.T ) + Chris_Ccc.Ptr * ( u_C.T 
             * u_C.R ) + Chris_Ccc.Pta * ( u_C.T * u_C.A ) + Chris_Ccc.Ptp 
             * ( u_C.T * u_C.P ) + Chris_Ccc.Prt * ( u_C.R * u_C.T ) 
             + Chris_Ccc.Prr * ( u_C.R * u_C.R ) + Chris_Ccc.Pra * ( u_C.R 
             * u_C.A ) + Chris_Ccc.Prp * ( u_C.R * u_C.P ) + Chris_Ccc.Pat 
             * ( u_C.A * u_C.T ) + Chris_Ccc.Par * ( u_C.A * u_C.R ) 
             + Chris_Ccc.Paa * ( u_C.A * u_C.A ) + Chris_Ccc.Pap * ( u_C.A 
             * u_C.P ) + Chris_Ccc.Ppt * ( u_C.P * u_C.T ) + Chris_Ccc.Ppr 
             * ( u_C.P * u_C.R ) + Chris_Ccc.Ppa * ( u_C.P * u_C.A ) 
             + Chris_Ccc.Ppp * ( u_C.P * u_C.P ) );

 u_C.T += ( accel_C.T * dtau );
 u_C.R += ( accel_C.R * dtau );
 u_C.A += ( accel_C.A * dtau );
 u_C.P += ( accel_C.P * dtau );

 r_C.T += ( u_C.T * dtau );
 r_C.R += ( u_C.R * dtau );
 r_C.A += ( u_C.A * dtau );
 r_C.P += ( u_C.P * dtau );

 t = r_C.T;
 r = r_C.R;
 a = r_C.A;
 p = r_C.P;

 ut = u_C.T;
 ur = u_C.R;
 ua = u_C.A;
 up = u_C.P;

/* Close GRPP Block */
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
