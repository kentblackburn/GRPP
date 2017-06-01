/* 
   	File:		example4.g
   	Author:		James Kent Blackburn 
	Date:		December 1993
   	Purpose:	Determine curvature for Kerr black hole

	Copyright (c) 1993-1998 
	By James Kent Blackburn
	All Rights Reserved
*/


/* Open  GRPP Block */
/* set up coordinates and indices used by GRPP */
/* COORDINATES: TtRrAaPp */
/* DIMENSION: 4 */

#include "example4_grpp.h" 

/* INDICES: IiJjKkLlMm */
/* TOTAL: 5 */

/* Close GRPP Block */

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

/* Open  GRPP Block */
 TENSOR_cc *g_cc ;
/* Close GRPP Block */
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


/* Open  GRPP Block */
/* assign all components of kerr metric tensor */
 (*g_cc).tt = gtt;

 (*g_cc).rr = grr;

 (*g_cc).aa = gaa;

 (*g_cc).pp = gpp;

 (*g_cc).pt = ( (*g_cc).tp = gtp );

 (*g_cc).rt = ( (*g_cc).tr = 0.0 );

 (*g_cc).at = ( (*g_cc).ta = 0.0 );

 (*g_cc).ar = ( (*g_cc).ra = 0.0 );

 (*g_cc).pr = ( (*g_cc).rp = 0.0 );

 (*g_cc).pa = ( (*g_cc).ap = 0.0 );

/* Close GRPP Block */
}

void kerr_christoffel(radius, theta, phi, Christoffel_Ccc)
double radius, theta, phi;

/* Open  GRPP Block */
 TENSOR_Ccc *Christoffel_Ccc ;
/* Close GRPP Block */
{

/* Open  GRPP Block */

 TENSOR_cc g_cc ;
 TENSOR_CC g_CC ;

 TENSOR_cc g_r_plus_eps_cc ;
 TENSOR_cc g_r_minus_eps_cc ;

 TENSOR_cc g_a_plus_eps_cc ;
 TENSOR_cc g_a_minus_eps_cc ;

 TENSOR_ccc _d_g_ccc ;

/* Close GRPP Block */
 void kerr_metric(), inverse();

 kerr_metric(radius, theta, phi, &g_cc);

/* Open  GRPP Block */
 array4x4[0][0] = g_cc.tt;
 array4x4[0][1] = g_cc.tr;
 array4x4[0][2] = g_cc.ta;
 array4x4[0][3] = g_cc.tp;
 array4x4[1][0] = g_cc.rt;
 array4x4[1][1] = g_cc.rr;
 array4x4[1][2] = g_cc.ra;
 array4x4[1][3] = g_cc.rp;
 array4x4[2][0] = g_cc.at;
 array4x4[2][1] = g_cc.ar;
 array4x4[2][2] = g_cc.aa;
 array4x4[2][3] = g_cc.ap;
 array4x4[3][0] = g_cc.pt;
 array4x4[3][1] = g_cc.pr;
 array4x4[3][2] = g_cc.pa;
 array4x4[3][3] = g_cc.pp;
/* Close GRPP Block */
 inverse(array4x4,4);

/* Open  GRPP Block */
 g_CC.TT = array4x4[0][0];
 g_CC.TR = array4x4[0][1];
 g_CC.TA = array4x4[0][2];
 g_CC.TP = array4x4[0][3];
 g_CC.RT = array4x4[1][0];
 g_CC.RR = array4x4[1][1];
 g_CC.RA = array4x4[1][2];
 g_CC.RP = array4x4[1][3];
 g_CC.AT = array4x4[2][0];
 g_CC.AR = array4x4[2][1];
 g_CC.AA = array4x4[2][2];
 g_CC.AP = array4x4[2][3];
 g_CC.PT = array4x4[3][0];
 g_CC.PR = array4x4[3][1];
 g_CC.PA = array4x4[3][2];
 g_CC.PP = array4x4[3][3];
/* Close GRPP Block */
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

/* Open  GRPP Block */

 _d_g_ccc.ttt = 0.0;
 _d_g_ccc.trt = 0.0;
 _d_g_ccc.tat = 0.0;
 _d_g_ccc.tpt = 0.0;
 _d_g_ccc.rtt = 0.0;
 _d_g_ccc.rrt = 0.0;
 _d_g_ccc.rat = 0.0;
 _d_g_ccc.rpt = 0.0;
 _d_g_ccc.att = 0.0;
 _d_g_ccc.art = 0.0;
 _d_g_ccc.aat = 0.0;
 _d_g_ccc.apt = 0.0;
 _d_g_ccc.ptt = 0.0;
 _d_g_ccc.prt = 0.0;
 _d_g_ccc.pat = 0.0;
 _d_g_ccc.ppt = 0.0;

 _d_g_ccc.ttr = ( g_r_plus_eps_cc.tt - g_r_minus_eps_cc.tt ) / ( 2 * EPS );
 _d_g_ccc.trr = ( g_r_plus_eps_cc.tr - g_r_minus_eps_cc.tr ) / ( 2 * EPS );
 _d_g_ccc.tar = ( g_r_plus_eps_cc.ta - g_r_minus_eps_cc.ta ) / ( 2 * EPS );
 _d_g_ccc.tpr = ( g_r_plus_eps_cc.tp - g_r_minus_eps_cc.tp ) / ( 2 * EPS );
 _d_g_ccc.rtr = ( g_r_plus_eps_cc.rt - g_r_minus_eps_cc.rt ) / ( 2 * EPS );
 _d_g_ccc.rrr = ( g_r_plus_eps_cc.rr - g_r_minus_eps_cc.rr ) / ( 2 * EPS );
 _d_g_ccc.rar = ( g_r_plus_eps_cc.ra - g_r_minus_eps_cc.ra ) / ( 2 * EPS );
 _d_g_ccc.rpr = ( g_r_plus_eps_cc.rp - g_r_minus_eps_cc.rp ) / ( 2 * EPS );
 _d_g_ccc.atr = ( g_r_plus_eps_cc.at - g_r_minus_eps_cc.at ) / ( 2 * EPS );
 _d_g_ccc.arr = ( g_r_plus_eps_cc.ar - g_r_minus_eps_cc.ar ) / ( 2 * EPS );
 _d_g_ccc.aar = ( g_r_plus_eps_cc.aa - g_r_minus_eps_cc.aa ) / ( 2 * EPS );
 _d_g_ccc.apr = ( g_r_plus_eps_cc.ap - g_r_minus_eps_cc.ap ) / ( 2 * EPS );
 _d_g_ccc.ptr = ( g_r_plus_eps_cc.pt - g_r_minus_eps_cc.pt ) / ( 2 * EPS );
 _d_g_ccc.prr = ( g_r_plus_eps_cc.pr - g_r_minus_eps_cc.pr ) / ( 2 * EPS );
 _d_g_ccc.par = ( g_r_plus_eps_cc.pa - g_r_minus_eps_cc.pa ) / ( 2 * EPS );
 _d_g_ccc.ppr = ( g_r_plus_eps_cc.pp - g_r_minus_eps_cc.pp ) / ( 2 * EPS );

 _d_g_ccc.tta = ( g_a_plus_eps_cc.tt - g_a_minus_eps_cc.tt ) / ( 2 * EPS );
 _d_g_ccc.tra = ( g_a_plus_eps_cc.tr - g_a_minus_eps_cc.tr ) / ( 2 * EPS );
 _d_g_ccc.taa = ( g_a_plus_eps_cc.ta - g_a_minus_eps_cc.ta ) / ( 2 * EPS );
 _d_g_ccc.tpa = ( g_a_plus_eps_cc.tp - g_a_minus_eps_cc.tp ) / ( 2 * EPS );
 _d_g_ccc.rta = ( g_a_plus_eps_cc.rt - g_a_minus_eps_cc.rt ) / ( 2 * EPS );
 _d_g_ccc.rra = ( g_a_plus_eps_cc.rr - g_a_minus_eps_cc.rr ) / ( 2 * EPS );
 _d_g_ccc.raa = ( g_a_plus_eps_cc.ra - g_a_minus_eps_cc.ra ) / ( 2 * EPS );
 _d_g_ccc.rpa = ( g_a_plus_eps_cc.rp - g_a_minus_eps_cc.rp ) / ( 2 * EPS );
 _d_g_ccc.ata = ( g_a_plus_eps_cc.at - g_a_minus_eps_cc.at ) / ( 2 * EPS );
 _d_g_ccc.ara = ( g_a_plus_eps_cc.ar - g_a_minus_eps_cc.ar ) / ( 2 * EPS );
 _d_g_ccc.aaa = ( g_a_plus_eps_cc.aa - g_a_minus_eps_cc.aa ) / ( 2 * EPS );
 _d_g_ccc.apa = ( g_a_plus_eps_cc.ap - g_a_minus_eps_cc.ap ) / ( 2 * EPS );
 _d_g_ccc.pta = ( g_a_plus_eps_cc.pt - g_a_minus_eps_cc.pt ) / ( 2 * EPS );
 _d_g_ccc.pra = ( g_a_plus_eps_cc.pr - g_a_minus_eps_cc.pr ) / ( 2 * EPS );
 _d_g_ccc.paa = ( g_a_plus_eps_cc.pa - g_a_minus_eps_cc.pa ) / ( 2 * EPS );
 _d_g_ccc.ppa = ( g_a_plus_eps_cc.pp - g_a_minus_eps_cc.pp ) / ( 2 * EPS );

 _d_g_ccc.ttp = 0.0;
 _d_g_ccc.trp = 0.0;
 _d_g_ccc.tap = 0.0;
 _d_g_ccc.tpp = 0.0;
 _d_g_ccc.rtp = 0.0;
 _d_g_ccc.rrp = 0.0;
 _d_g_ccc.rap = 0.0;
 _d_g_ccc.rpp = 0.0;
 _d_g_ccc.atp = 0.0;
 _d_g_ccc.arp = 0.0;
 _d_g_ccc.aap = 0.0;
 _d_g_ccc.app = 0.0;
 _d_g_ccc.ptp = 0.0;
 _d_g_ccc.prp = 0.0;
 _d_g_ccc.pap = 0.0;
 _d_g_ccc.ppp = 0.0;

 (*Christoffel_Ccc).Ttt = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                          - _d_g_ccc.ttt ) + g_CC.TR * ( _d_g_ccc.rtt 
                          + _d_g_ccc.rtt - _d_g_ccc.ttr ) + g_CC.TA 
                          * ( _d_g_ccc.att + _d_g_ccc.att - _d_g_ccc.tta ) 
                          + g_CC.TP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                          - _d_g_ccc.ttp ) ) );
 (*Christoffel_Ccc).Ttr = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                          - _d_g_ccc.rtt ) + g_CC.TR * ( _d_g_ccc.rtr 
                          + _d_g_ccc.rrt - _d_g_ccc.rtr ) + g_CC.TA 
                          * ( _d_g_ccc.atr + _d_g_ccc.art - _d_g_ccc.rta ) 
                          + g_CC.TP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                          - _d_g_ccc.rtp ) ) );
 (*Christoffel_Ccc).Tta = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                          - _d_g_ccc.att ) + g_CC.TR * ( _d_g_ccc.rta 
                          + _d_g_ccc.rat - _d_g_ccc.atr ) + g_CC.TA 
                          * ( _d_g_ccc.ata + _d_g_ccc.aat - _d_g_ccc.ata ) 
                          + g_CC.TP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                          - _d_g_ccc.atp ) ) );
 (*Christoffel_Ccc).Ttp = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                          - _d_g_ccc.ptt ) + g_CC.TR * ( _d_g_ccc.rtp 
                          + _d_g_ccc.rpt - _d_g_ccc.ptr ) + g_CC.TA 
                          * ( _d_g_ccc.atp + _d_g_ccc.apt - _d_g_ccc.pta ) 
                          + g_CC.TP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                          - _d_g_ccc.ptp ) ) );
 (*Christoffel_Ccc).Trt = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                          - _d_g_ccc.trt ) + g_CC.TR * ( _d_g_ccc.rrt 
                          + _d_g_ccc.rtr - _d_g_ccc.trr ) + g_CC.TA 
                          * ( _d_g_ccc.art + _d_g_ccc.atr - _d_g_ccc.tra ) 
                          + g_CC.TP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                          - _d_g_ccc.trp ) ) );
 (*Christoffel_Ccc).Trr = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                          - _d_g_ccc.rrt ) + g_CC.TR * ( _d_g_ccc.rrr 
                          + _d_g_ccc.rrr - _d_g_ccc.rrr ) + g_CC.TA 
                          * ( _d_g_ccc.arr + _d_g_ccc.arr - _d_g_ccc.rra ) 
                          + g_CC.TP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                          - _d_g_ccc.rrp ) ) );
 (*Christoffel_Ccc).Tra = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                          - _d_g_ccc.art ) + g_CC.TR * ( _d_g_ccc.rra 
                          + _d_g_ccc.rar - _d_g_ccc.arr ) + g_CC.TA 
                          * ( _d_g_ccc.ara + _d_g_ccc.aar - _d_g_ccc.ara ) 
                          + g_CC.TP * ( _d_g_ccc.pra + _d_g_ccc.par 
                          - _d_g_ccc.arp ) ) );
 (*Christoffel_Ccc).Trp = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                          - _d_g_ccc.prt ) + g_CC.TR * ( _d_g_ccc.rrp 
                          + _d_g_ccc.rpr - _d_g_ccc.prr ) + g_CC.TA 
                          * ( _d_g_ccc.arp + _d_g_ccc.apr - _d_g_ccc.pra ) 
                          + g_CC.TP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                          - _d_g_ccc.prp ) ) );
 (*Christoffel_Ccc).Tat = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                          - _d_g_ccc.tat ) + g_CC.TR * ( _d_g_ccc.rat 
                          + _d_g_ccc.rta - _d_g_ccc.tar ) + g_CC.TA 
                          * ( _d_g_ccc.aat + _d_g_ccc.ata - _d_g_ccc.taa ) 
                          + g_CC.TP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                          - _d_g_ccc.tap ) ) );
 (*Christoffel_Ccc).Tar = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                          - _d_g_ccc.rat ) + g_CC.TR * ( _d_g_ccc.rar 
                          + _d_g_ccc.rra - _d_g_ccc.rar ) + g_CC.TA 
                          * ( _d_g_ccc.aar + _d_g_ccc.ara - _d_g_ccc.raa ) 
                          + g_CC.TP * ( _d_g_ccc.par + _d_g_ccc.pra 
                          - _d_g_ccc.rap ) ) );
 (*Christoffel_Ccc).Taa = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                          - _d_g_ccc.aat ) + g_CC.TR * ( _d_g_ccc.raa 
                          + _d_g_ccc.raa - _d_g_ccc.aar ) + g_CC.TA 
                          * ( _d_g_ccc.aaa + _d_g_ccc.aaa - _d_g_ccc.aaa ) 
                          + g_CC.TP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                          - _d_g_ccc.aap ) ) );
 (*Christoffel_Ccc).Tap = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                          - _d_g_ccc.pat ) + g_CC.TR * ( _d_g_ccc.rap 
                          + _d_g_ccc.rpa - _d_g_ccc.par ) + g_CC.TA 
                          * ( _d_g_ccc.aap + _d_g_ccc.apa - _d_g_ccc.paa ) 
                          + g_CC.TP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                          - _d_g_ccc.pap ) ) );
 (*Christoffel_Ccc).Tpt = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                          - _d_g_ccc.tpt ) + g_CC.TR * ( _d_g_ccc.rpt 
                          + _d_g_ccc.rtp - _d_g_ccc.tpr ) + g_CC.TA 
                          * ( _d_g_ccc.apt + _d_g_ccc.atp - _d_g_ccc.tpa ) 
                          + g_CC.TP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                          - _d_g_ccc.tpp ) ) );
 (*Christoffel_Ccc).Tpr = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                          - _d_g_ccc.rpt ) + g_CC.TR * ( _d_g_ccc.rpr 
                          + _d_g_ccc.rrp - _d_g_ccc.rpr ) + g_CC.TA 
                          * ( _d_g_ccc.apr + _d_g_ccc.arp - _d_g_ccc.rpa ) 
                          + g_CC.TP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                          - _d_g_ccc.rpp ) ) );
 (*Christoffel_Ccc).Tpa = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                          - _d_g_ccc.apt ) + g_CC.TR * ( _d_g_ccc.rpa 
                          + _d_g_ccc.rap - _d_g_ccc.apr ) + g_CC.TA 
                          * ( _d_g_ccc.apa + _d_g_ccc.aap - _d_g_ccc.apa ) 
                          + g_CC.TP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                          - _d_g_ccc.app ) ) );
 (*Christoffel_Ccc).Tpp = ( 0.5 * ( g_CC.TT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                          - _d_g_ccc.ppt ) + g_CC.TR * ( _d_g_ccc.rpp 
                          + _d_g_ccc.rpp - _d_g_ccc.ppr ) + g_CC.TA 
                          * ( _d_g_ccc.app + _d_g_ccc.app - _d_g_ccc.ppa ) 
                          + g_CC.TP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                          - _d_g_ccc.ppp ) ) );
 (*Christoffel_Ccc).Rtt = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                          - _d_g_ccc.ttt ) + g_CC.RR * ( _d_g_ccc.rtt 
                          + _d_g_ccc.rtt - _d_g_ccc.ttr ) + g_CC.RA 
                          * ( _d_g_ccc.att + _d_g_ccc.att - _d_g_ccc.tta ) 
                          + g_CC.RP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                          - _d_g_ccc.ttp ) ) );
 (*Christoffel_Ccc).Rtr = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                          - _d_g_ccc.rtt ) + g_CC.RR * ( _d_g_ccc.rtr 
                          + _d_g_ccc.rrt - _d_g_ccc.rtr ) + g_CC.RA 
                          * ( _d_g_ccc.atr + _d_g_ccc.art - _d_g_ccc.rta ) 
                          + g_CC.RP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                          - _d_g_ccc.rtp ) ) );
 (*Christoffel_Ccc).Rta = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                          - _d_g_ccc.att ) + g_CC.RR * ( _d_g_ccc.rta 
                          + _d_g_ccc.rat - _d_g_ccc.atr ) + g_CC.RA 
                          * ( _d_g_ccc.ata + _d_g_ccc.aat - _d_g_ccc.ata ) 
                          + g_CC.RP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                          - _d_g_ccc.atp ) ) );
 (*Christoffel_Ccc).Rtp = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                          - _d_g_ccc.ptt ) + g_CC.RR * ( _d_g_ccc.rtp 
                          + _d_g_ccc.rpt - _d_g_ccc.ptr ) + g_CC.RA 
                          * ( _d_g_ccc.atp + _d_g_ccc.apt - _d_g_ccc.pta ) 
                          + g_CC.RP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                          - _d_g_ccc.ptp ) ) );
 (*Christoffel_Ccc).Rrt = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                          - _d_g_ccc.trt ) + g_CC.RR * ( _d_g_ccc.rrt 
                          + _d_g_ccc.rtr - _d_g_ccc.trr ) + g_CC.RA 
                          * ( _d_g_ccc.art + _d_g_ccc.atr - _d_g_ccc.tra ) 
                          + g_CC.RP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                          - _d_g_ccc.trp ) ) );
 (*Christoffel_Ccc).Rrr = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                          - _d_g_ccc.rrt ) + g_CC.RR * ( _d_g_ccc.rrr 
                          + _d_g_ccc.rrr - _d_g_ccc.rrr ) + g_CC.RA 
                          * ( _d_g_ccc.arr + _d_g_ccc.arr - _d_g_ccc.rra ) 
                          + g_CC.RP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                          - _d_g_ccc.rrp ) ) );
 (*Christoffel_Ccc).Rra = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                          - _d_g_ccc.art ) + g_CC.RR * ( _d_g_ccc.rra 
                          + _d_g_ccc.rar - _d_g_ccc.arr ) + g_CC.RA 
                          * ( _d_g_ccc.ara + _d_g_ccc.aar - _d_g_ccc.ara ) 
                          + g_CC.RP * ( _d_g_ccc.pra + _d_g_ccc.par 
                          - _d_g_ccc.arp ) ) );
 (*Christoffel_Ccc).Rrp = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                          - _d_g_ccc.prt ) + g_CC.RR * ( _d_g_ccc.rrp 
                          + _d_g_ccc.rpr - _d_g_ccc.prr ) + g_CC.RA 
                          * ( _d_g_ccc.arp + _d_g_ccc.apr - _d_g_ccc.pra ) 
                          + g_CC.RP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                          - _d_g_ccc.prp ) ) );
 (*Christoffel_Ccc).Rat = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                          - _d_g_ccc.tat ) + g_CC.RR * ( _d_g_ccc.rat 
                          + _d_g_ccc.rta - _d_g_ccc.tar ) + g_CC.RA 
                          * ( _d_g_ccc.aat + _d_g_ccc.ata - _d_g_ccc.taa ) 
                          + g_CC.RP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                          - _d_g_ccc.tap ) ) );
 (*Christoffel_Ccc).Rar = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                          - _d_g_ccc.rat ) + g_CC.RR * ( _d_g_ccc.rar 
                          + _d_g_ccc.rra - _d_g_ccc.rar ) + g_CC.RA 
                          * ( _d_g_ccc.aar + _d_g_ccc.ara - _d_g_ccc.raa ) 
                          + g_CC.RP * ( _d_g_ccc.par + _d_g_ccc.pra 
                          - _d_g_ccc.rap ) ) );
 (*Christoffel_Ccc).Raa = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                          - _d_g_ccc.aat ) + g_CC.RR * ( _d_g_ccc.raa 
                          + _d_g_ccc.raa - _d_g_ccc.aar ) + g_CC.RA 
                          * ( _d_g_ccc.aaa + _d_g_ccc.aaa - _d_g_ccc.aaa ) 
                          + g_CC.RP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                          - _d_g_ccc.aap ) ) );
 (*Christoffel_Ccc).Rap = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                          - _d_g_ccc.pat ) + g_CC.RR * ( _d_g_ccc.rap 
                          + _d_g_ccc.rpa - _d_g_ccc.par ) + g_CC.RA 
                          * ( _d_g_ccc.aap + _d_g_ccc.apa - _d_g_ccc.paa ) 
                          + g_CC.RP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                          - _d_g_ccc.pap ) ) );
 (*Christoffel_Ccc).Rpt = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                          - _d_g_ccc.tpt ) + g_CC.RR * ( _d_g_ccc.rpt 
                          + _d_g_ccc.rtp - _d_g_ccc.tpr ) + g_CC.RA 
                          * ( _d_g_ccc.apt + _d_g_ccc.atp - _d_g_ccc.tpa ) 
                          + g_CC.RP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                          - _d_g_ccc.tpp ) ) );
 (*Christoffel_Ccc).Rpr = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                          - _d_g_ccc.rpt ) + g_CC.RR * ( _d_g_ccc.rpr 
                          + _d_g_ccc.rrp - _d_g_ccc.rpr ) + g_CC.RA 
                          * ( _d_g_ccc.apr + _d_g_ccc.arp - _d_g_ccc.rpa ) 
                          + g_CC.RP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                          - _d_g_ccc.rpp ) ) );
 (*Christoffel_Ccc).Rpa = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                          - _d_g_ccc.apt ) + g_CC.RR * ( _d_g_ccc.rpa 
                          + _d_g_ccc.rap - _d_g_ccc.apr ) + g_CC.RA 
                          * ( _d_g_ccc.apa + _d_g_ccc.aap - _d_g_ccc.apa ) 
                          + g_CC.RP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                          - _d_g_ccc.app ) ) );
 (*Christoffel_Ccc).Rpp = ( 0.5 * ( g_CC.RT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                          - _d_g_ccc.ppt ) + g_CC.RR * ( _d_g_ccc.rpp 
                          + _d_g_ccc.rpp - _d_g_ccc.ppr ) + g_CC.RA 
                          * ( _d_g_ccc.app + _d_g_ccc.app - _d_g_ccc.ppa ) 
                          + g_CC.RP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                          - _d_g_ccc.ppp ) ) );
 (*Christoffel_Ccc).Att = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                          - _d_g_ccc.ttt ) + g_CC.AR * ( _d_g_ccc.rtt 
                          + _d_g_ccc.rtt - _d_g_ccc.ttr ) + g_CC.AA 
                          * ( _d_g_ccc.att + _d_g_ccc.att - _d_g_ccc.tta ) 
                          + g_CC.AP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                          - _d_g_ccc.ttp ) ) );
 (*Christoffel_Ccc).Atr = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                          - _d_g_ccc.rtt ) + g_CC.AR * ( _d_g_ccc.rtr 
                          + _d_g_ccc.rrt - _d_g_ccc.rtr ) + g_CC.AA 
                          * ( _d_g_ccc.atr + _d_g_ccc.art - _d_g_ccc.rta ) 
                          + g_CC.AP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                          - _d_g_ccc.rtp ) ) );
 (*Christoffel_Ccc).Ata = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                          - _d_g_ccc.att ) + g_CC.AR * ( _d_g_ccc.rta 
                          + _d_g_ccc.rat - _d_g_ccc.atr ) + g_CC.AA 
                          * ( _d_g_ccc.ata + _d_g_ccc.aat - _d_g_ccc.ata ) 
                          + g_CC.AP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                          - _d_g_ccc.atp ) ) );
 (*Christoffel_Ccc).Atp = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                          - _d_g_ccc.ptt ) + g_CC.AR * ( _d_g_ccc.rtp 
                          + _d_g_ccc.rpt - _d_g_ccc.ptr ) + g_CC.AA 
                          * ( _d_g_ccc.atp + _d_g_ccc.apt - _d_g_ccc.pta ) 
                          + g_CC.AP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                          - _d_g_ccc.ptp ) ) );
 (*Christoffel_Ccc).Art = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                          - _d_g_ccc.trt ) + g_CC.AR * ( _d_g_ccc.rrt 
                          + _d_g_ccc.rtr - _d_g_ccc.trr ) + g_CC.AA 
                          * ( _d_g_ccc.art + _d_g_ccc.atr - _d_g_ccc.tra ) 
                          + g_CC.AP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                          - _d_g_ccc.trp ) ) );
 (*Christoffel_Ccc).Arr = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                          - _d_g_ccc.rrt ) + g_CC.AR * ( _d_g_ccc.rrr 
                          + _d_g_ccc.rrr - _d_g_ccc.rrr ) + g_CC.AA 
                          * ( _d_g_ccc.arr + _d_g_ccc.arr - _d_g_ccc.rra ) 
                          + g_CC.AP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                          - _d_g_ccc.rrp ) ) );
 (*Christoffel_Ccc).Ara = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                          - _d_g_ccc.art ) + g_CC.AR * ( _d_g_ccc.rra 
                          + _d_g_ccc.rar - _d_g_ccc.arr ) + g_CC.AA 
                          * ( _d_g_ccc.ara + _d_g_ccc.aar - _d_g_ccc.ara ) 
                          + g_CC.AP * ( _d_g_ccc.pra + _d_g_ccc.par 
                          - _d_g_ccc.arp ) ) );
 (*Christoffel_Ccc).Arp = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                          - _d_g_ccc.prt ) + g_CC.AR * ( _d_g_ccc.rrp 
                          + _d_g_ccc.rpr - _d_g_ccc.prr ) + g_CC.AA 
                          * ( _d_g_ccc.arp + _d_g_ccc.apr - _d_g_ccc.pra ) 
                          + g_CC.AP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                          - _d_g_ccc.prp ) ) );
 (*Christoffel_Ccc).Aat = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                          - _d_g_ccc.tat ) + g_CC.AR * ( _d_g_ccc.rat 
                          + _d_g_ccc.rta - _d_g_ccc.tar ) + g_CC.AA 
                          * ( _d_g_ccc.aat + _d_g_ccc.ata - _d_g_ccc.taa ) 
                          + g_CC.AP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                          - _d_g_ccc.tap ) ) );
 (*Christoffel_Ccc).Aar = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                          - _d_g_ccc.rat ) + g_CC.AR * ( _d_g_ccc.rar 
                          + _d_g_ccc.rra - _d_g_ccc.rar ) + g_CC.AA 
                          * ( _d_g_ccc.aar + _d_g_ccc.ara - _d_g_ccc.raa ) 
                          + g_CC.AP * ( _d_g_ccc.par + _d_g_ccc.pra 
                          - _d_g_ccc.rap ) ) );
 (*Christoffel_Ccc).Aaa = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                          - _d_g_ccc.aat ) + g_CC.AR * ( _d_g_ccc.raa 
                          + _d_g_ccc.raa - _d_g_ccc.aar ) + g_CC.AA 
                          * ( _d_g_ccc.aaa + _d_g_ccc.aaa - _d_g_ccc.aaa ) 
                          + g_CC.AP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                          - _d_g_ccc.aap ) ) );
 (*Christoffel_Ccc).Aap = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                          - _d_g_ccc.pat ) + g_CC.AR * ( _d_g_ccc.rap 
                          + _d_g_ccc.rpa - _d_g_ccc.par ) + g_CC.AA 
                          * ( _d_g_ccc.aap + _d_g_ccc.apa - _d_g_ccc.paa ) 
                          + g_CC.AP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                          - _d_g_ccc.pap ) ) );
 (*Christoffel_Ccc).Apt = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                          - _d_g_ccc.tpt ) + g_CC.AR * ( _d_g_ccc.rpt 
                          + _d_g_ccc.rtp - _d_g_ccc.tpr ) + g_CC.AA 
                          * ( _d_g_ccc.apt + _d_g_ccc.atp - _d_g_ccc.tpa ) 
                          + g_CC.AP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                          - _d_g_ccc.tpp ) ) );
 (*Christoffel_Ccc).Apr = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                          - _d_g_ccc.rpt ) + g_CC.AR * ( _d_g_ccc.rpr 
                          + _d_g_ccc.rrp - _d_g_ccc.rpr ) + g_CC.AA 
                          * ( _d_g_ccc.apr + _d_g_ccc.arp - _d_g_ccc.rpa ) 
                          + g_CC.AP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                          - _d_g_ccc.rpp ) ) );
 (*Christoffel_Ccc).Apa = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                          - _d_g_ccc.apt ) + g_CC.AR * ( _d_g_ccc.rpa 
                          + _d_g_ccc.rap - _d_g_ccc.apr ) + g_CC.AA 
                          * ( _d_g_ccc.apa + _d_g_ccc.aap - _d_g_ccc.apa ) 
                          + g_CC.AP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                          - _d_g_ccc.app ) ) );
 (*Christoffel_Ccc).App = ( 0.5 * ( g_CC.AT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                          - _d_g_ccc.ppt ) + g_CC.AR * ( _d_g_ccc.rpp 
                          + _d_g_ccc.rpp - _d_g_ccc.ppr ) + g_CC.AA 
                          * ( _d_g_ccc.app + _d_g_ccc.app - _d_g_ccc.ppa ) 
                          + g_CC.AP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                          - _d_g_ccc.ppp ) ) );
 (*Christoffel_Ccc).Ptt = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.ttt + _d_g_ccc.ttt 
                          - _d_g_ccc.ttt ) + g_CC.PR * ( _d_g_ccc.rtt 
                          + _d_g_ccc.rtt - _d_g_ccc.ttr ) + g_CC.PA 
                          * ( _d_g_ccc.att + _d_g_ccc.att - _d_g_ccc.tta ) 
                          + g_CC.PP * ( _d_g_ccc.ptt + _d_g_ccc.ptt 
                          - _d_g_ccc.ttp ) ) );
 (*Christoffel_Ccc).Ptr = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.ttr + _d_g_ccc.trt 
                          - _d_g_ccc.rtt ) + g_CC.PR * ( _d_g_ccc.rtr 
                          + _d_g_ccc.rrt - _d_g_ccc.rtr ) + g_CC.PA 
                          * ( _d_g_ccc.atr + _d_g_ccc.art - _d_g_ccc.rta ) 
                          + g_CC.PP * ( _d_g_ccc.ptr + _d_g_ccc.prt 
                          - _d_g_ccc.rtp ) ) );
 (*Christoffel_Ccc).Pta = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tta + _d_g_ccc.tat 
                          - _d_g_ccc.att ) + g_CC.PR * ( _d_g_ccc.rta 
                          + _d_g_ccc.rat - _d_g_ccc.atr ) + g_CC.PA 
                          * ( _d_g_ccc.ata + _d_g_ccc.aat - _d_g_ccc.ata ) 
                          + g_CC.PP * ( _d_g_ccc.pta + _d_g_ccc.pat 
                          - _d_g_ccc.atp ) ) );
 (*Christoffel_Ccc).Ptp = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.ttp + _d_g_ccc.tpt 
                          - _d_g_ccc.ptt ) + g_CC.PR * ( _d_g_ccc.rtp 
                          + _d_g_ccc.rpt - _d_g_ccc.ptr ) + g_CC.PA 
                          * ( _d_g_ccc.atp + _d_g_ccc.apt - _d_g_ccc.pta ) 
                          + g_CC.PP * ( _d_g_ccc.ptp + _d_g_ccc.ppt 
                          - _d_g_ccc.ptp ) ) );
 (*Christoffel_Ccc).Prt = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.trt + _d_g_ccc.ttr 
                          - _d_g_ccc.trt ) + g_CC.PR * ( _d_g_ccc.rrt 
                          + _d_g_ccc.rtr - _d_g_ccc.trr ) + g_CC.PA 
                          * ( _d_g_ccc.art + _d_g_ccc.atr - _d_g_ccc.tra ) 
                          + g_CC.PP * ( _d_g_ccc.prt + _d_g_ccc.ptr 
                          - _d_g_ccc.trp ) ) );
 (*Christoffel_Ccc).Prr = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.trr + _d_g_ccc.trr 
                          - _d_g_ccc.rrt ) + g_CC.PR * ( _d_g_ccc.rrr 
                          + _d_g_ccc.rrr - _d_g_ccc.rrr ) + g_CC.PA 
                          * ( _d_g_ccc.arr + _d_g_ccc.arr - _d_g_ccc.rra ) 
                          + g_CC.PP * ( _d_g_ccc.prr + _d_g_ccc.prr 
                          - _d_g_ccc.rrp ) ) );
 (*Christoffel_Ccc).Pra = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tra + _d_g_ccc.tar 
                          - _d_g_ccc.art ) + g_CC.PR * ( _d_g_ccc.rra 
                          + _d_g_ccc.rar - _d_g_ccc.arr ) + g_CC.PA 
                          * ( _d_g_ccc.ara + _d_g_ccc.aar - _d_g_ccc.ara ) 
                          + g_CC.PP * ( _d_g_ccc.pra + _d_g_ccc.par 
                          - _d_g_ccc.arp ) ) );
 (*Christoffel_Ccc).Prp = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.trp + _d_g_ccc.tpr 
                          - _d_g_ccc.prt ) + g_CC.PR * ( _d_g_ccc.rrp 
                          + _d_g_ccc.rpr - _d_g_ccc.prr ) + g_CC.PA 
                          * ( _d_g_ccc.arp + _d_g_ccc.apr - _d_g_ccc.pra ) 
                          + g_CC.PP * ( _d_g_ccc.prp + _d_g_ccc.ppr 
                          - _d_g_ccc.prp ) ) );
 (*Christoffel_Ccc).Pat = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tat + _d_g_ccc.tta 
                          - _d_g_ccc.tat ) + g_CC.PR * ( _d_g_ccc.rat 
                          + _d_g_ccc.rta - _d_g_ccc.tar ) + g_CC.PA 
                          * ( _d_g_ccc.aat + _d_g_ccc.ata - _d_g_ccc.taa ) 
                          + g_CC.PP * ( _d_g_ccc.pat + _d_g_ccc.pta 
                          - _d_g_ccc.tap ) ) );
 (*Christoffel_Ccc).Par = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tar + _d_g_ccc.tra 
                          - _d_g_ccc.rat ) + g_CC.PR * ( _d_g_ccc.rar 
                          + _d_g_ccc.rra - _d_g_ccc.rar ) + g_CC.PA 
                          * ( _d_g_ccc.aar + _d_g_ccc.ara - _d_g_ccc.raa ) 
                          + g_CC.PP * ( _d_g_ccc.par + _d_g_ccc.pra 
                          - _d_g_ccc.rap ) ) );
 (*Christoffel_Ccc).Paa = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.taa + _d_g_ccc.taa 
                          - _d_g_ccc.aat ) + g_CC.PR * ( _d_g_ccc.raa 
                          + _d_g_ccc.raa - _d_g_ccc.aar ) + g_CC.PA 
                          * ( _d_g_ccc.aaa + _d_g_ccc.aaa - _d_g_ccc.aaa ) 
                          + g_CC.PP * ( _d_g_ccc.paa + _d_g_ccc.paa 
                          - _d_g_ccc.aap ) ) );
 (*Christoffel_Ccc).Pap = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tap + _d_g_ccc.tpa 
                          - _d_g_ccc.pat ) + g_CC.PR * ( _d_g_ccc.rap 
                          + _d_g_ccc.rpa - _d_g_ccc.par ) + g_CC.PA 
                          * ( _d_g_ccc.aap + _d_g_ccc.apa - _d_g_ccc.paa ) 
                          + g_CC.PP * ( _d_g_ccc.pap + _d_g_ccc.ppa 
                          - _d_g_ccc.pap ) ) );
 (*Christoffel_Ccc).Ppt = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpt + _d_g_ccc.ttp 
                          - _d_g_ccc.tpt ) + g_CC.PR * ( _d_g_ccc.rpt 
                          + _d_g_ccc.rtp - _d_g_ccc.tpr ) + g_CC.PA 
                          * ( _d_g_ccc.apt + _d_g_ccc.atp - _d_g_ccc.tpa ) 
                          + g_CC.PP * ( _d_g_ccc.ppt + _d_g_ccc.ptp 
                          - _d_g_ccc.tpp ) ) );
 (*Christoffel_Ccc).Ppr = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpr + _d_g_ccc.trp 
                          - _d_g_ccc.rpt ) + g_CC.PR * ( _d_g_ccc.rpr 
                          + _d_g_ccc.rrp - _d_g_ccc.rpr ) + g_CC.PA 
                          * ( _d_g_ccc.apr + _d_g_ccc.arp - _d_g_ccc.rpa ) 
                          + g_CC.PP * ( _d_g_ccc.ppr + _d_g_ccc.prp 
                          - _d_g_ccc.rpp ) ) );
 (*Christoffel_Ccc).Ppa = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpa + _d_g_ccc.tap 
                          - _d_g_ccc.apt ) + g_CC.PR * ( _d_g_ccc.rpa 
                          + _d_g_ccc.rap - _d_g_ccc.apr ) + g_CC.PA 
                          * ( _d_g_ccc.apa + _d_g_ccc.aap - _d_g_ccc.apa ) 
                          + g_CC.PP * ( _d_g_ccc.ppa + _d_g_ccc.pap 
                          - _d_g_ccc.app ) ) );
 (*Christoffel_Ccc).Ppp = ( 0.5 * ( g_CC.PT * ( _d_g_ccc.tpp + _d_g_ccc.tpp 
                          - _d_g_ccc.ppt ) + g_CC.PR * ( _d_g_ccc.rpp 
                          + _d_g_ccc.rpp - _d_g_ccc.ppr ) + g_CC.PA 
                          * ( _d_g_ccc.app + _d_g_ccc.app - _d_g_ccc.ppa ) 
                          + g_CC.PP * ( _d_g_ccc.ppp + _d_g_ccc.ppp 
                          - _d_g_ccc.ppp ) ) );

/* Close GRPP Block */
}

double kerr_curvature(radius, theta, phi)
double radius, theta, phi;
{
 double R;

/* Open  GRPP Block */

 TENSOR_cc g_cc ;
 TENSOR_CC g_CC ;
 TENSOR_cc R_cc ;

 TENSOR_Ccc Chr_Ccc ;

 TENSOR_Ccc Chr_r_plus_eps_Ccc ;
 TENSOR_Ccc Chr_r_minus_eps_Ccc ;

 TENSOR_Ccc Chr_a_plus_eps_Ccc ;
 TENSOR_Ccc Chr_a_minus_eps_Ccc ;

 TENSOR_Cccc R_Cccc ;
 TENSOR_Cccc _d_Chr_Cccc ;

/* Close GRPP Block */
 void kerr_metric(), kerr_christoffel(), inverse();

 kerr_metric(radius, theta, phi, &g_cc);

/* Open  GRPP Block */
 array4x4[0][0] = g_cc.tt;
 array4x4[0][1] = g_cc.tr;
 array4x4[0][2] = g_cc.ta;
 array4x4[0][3] = g_cc.tp;
 array4x4[1][0] = g_cc.rt;
 array4x4[1][1] = g_cc.rr;
 array4x4[1][2] = g_cc.ra;
 array4x4[1][3] = g_cc.rp;
 array4x4[2][0] = g_cc.at;
 array4x4[2][1] = g_cc.ar;
 array4x4[2][2] = g_cc.aa;
 array4x4[2][3] = g_cc.ap;
 array4x4[3][0] = g_cc.pt;
 array4x4[3][1] = g_cc.pr;
 array4x4[3][2] = g_cc.pa;
 array4x4[3][3] = g_cc.pp;
/* Close GRPP Block */
 inverse(array4x4,4);

/* Open  GRPP Block */
 g_CC.TT = array4x4[0][0];
 g_CC.TR = array4x4[0][1];
 g_CC.TA = array4x4[0][2];
 g_CC.TP = array4x4[0][3];
 g_CC.RT = array4x4[1][0];
 g_CC.RR = array4x4[1][1];
 g_CC.RA = array4x4[1][2];
 g_CC.RP = array4x4[1][3];
 g_CC.AT = array4x4[2][0];
 g_CC.AR = array4x4[2][1];
 g_CC.AA = array4x4[2][2];
 g_CC.AP = array4x4[2][3];
 g_CC.PT = array4x4[3][0];
 g_CC.PR = array4x4[3][1];
 g_CC.PA = array4x4[3][2];
 g_CC.PP = array4x4[3][3];
/* Close GRPP Block */
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

/* Open  GRPP Block */

 _d_Chr_Cccc.Tttt = 0.0;
 _d_Chr_Cccc.Ttrt = 0.0;
 _d_Chr_Cccc.Ttat = 0.0;
 _d_Chr_Cccc.Ttpt = 0.0;
 _d_Chr_Cccc.Trtt = 0.0;
 _d_Chr_Cccc.Trrt = 0.0;
 _d_Chr_Cccc.Trat = 0.0;
 _d_Chr_Cccc.Trpt = 0.0;
 _d_Chr_Cccc.Tatt = 0.0;
 _d_Chr_Cccc.Tart = 0.0;
 _d_Chr_Cccc.Taat = 0.0;
 _d_Chr_Cccc.Tapt = 0.0;
 _d_Chr_Cccc.Tptt = 0.0;
 _d_Chr_Cccc.Tprt = 0.0;
 _d_Chr_Cccc.Tpat = 0.0;
 _d_Chr_Cccc.Tppt = 0.0;
 _d_Chr_Cccc.Rttt = 0.0;
 _d_Chr_Cccc.Rtrt = 0.0;
 _d_Chr_Cccc.Rtat = 0.0;
 _d_Chr_Cccc.Rtpt = 0.0;
 _d_Chr_Cccc.Rrtt = 0.0;
 _d_Chr_Cccc.Rrrt = 0.0;
 _d_Chr_Cccc.Rrat = 0.0;
 _d_Chr_Cccc.Rrpt = 0.0;
 _d_Chr_Cccc.Ratt = 0.0;
 _d_Chr_Cccc.Rart = 0.0;
 _d_Chr_Cccc.Raat = 0.0;
 _d_Chr_Cccc.Rapt = 0.0;
 _d_Chr_Cccc.Rptt = 0.0;
 _d_Chr_Cccc.Rprt = 0.0;
 _d_Chr_Cccc.Rpat = 0.0;
 _d_Chr_Cccc.Rppt = 0.0;
 _d_Chr_Cccc.Attt = 0.0;
 _d_Chr_Cccc.Atrt = 0.0;
 _d_Chr_Cccc.Atat = 0.0;
 _d_Chr_Cccc.Atpt = 0.0;
 _d_Chr_Cccc.Artt = 0.0;
 _d_Chr_Cccc.Arrt = 0.0;
 _d_Chr_Cccc.Arat = 0.0;
 _d_Chr_Cccc.Arpt = 0.0;
 _d_Chr_Cccc.Aatt = 0.0;
 _d_Chr_Cccc.Aart = 0.0;
 _d_Chr_Cccc.Aaat = 0.0;
 _d_Chr_Cccc.Aapt = 0.0;
 _d_Chr_Cccc.Aptt = 0.0;
 _d_Chr_Cccc.Aprt = 0.0;
 _d_Chr_Cccc.Apat = 0.0;
 _d_Chr_Cccc.Appt = 0.0;
 _d_Chr_Cccc.Pttt = 0.0;
 _d_Chr_Cccc.Ptrt = 0.0;
 _d_Chr_Cccc.Ptat = 0.0;
 _d_Chr_Cccc.Ptpt = 0.0;
 _d_Chr_Cccc.Prtt = 0.0;
 _d_Chr_Cccc.Prrt = 0.0;
 _d_Chr_Cccc.Prat = 0.0;
 _d_Chr_Cccc.Prpt = 0.0;
 _d_Chr_Cccc.Patt = 0.0;
 _d_Chr_Cccc.Part = 0.0;
 _d_Chr_Cccc.Paat = 0.0;
 _d_Chr_Cccc.Papt = 0.0;
 _d_Chr_Cccc.Pptt = 0.0;
 _d_Chr_Cccc.Pprt = 0.0;
 _d_Chr_Cccc.Ppat = 0.0;
 _d_Chr_Cccc.Pppt = 0.0;

 _d_Chr_Cccc.Tttr = ( Chr_r_plus_eps_Ccc.Ttt - Chr_r_minus_eps_Ccc.Ttt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ttrr = ( Chr_r_plus_eps_Ccc.Ttr - Chr_r_minus_eps_Ccc.Ttr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ttar = ( Chr_r_plus_eps_Ccc.Tta - Chr_r_minus_eps_Ccc.Tta ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ttpr = ( Chr_r_plus_eps_Ccc.Ttp - Chr_r_minus_eps_Ccc.Ttp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trtr = ( Chr_r_plus_eps_Ccc.Trt - Chr_r_minus_eps_Ccc.Trt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trrr = ( Chr_r_plus_eps_Ccc.Trr - Chr_r_minus_eps_Ccc.Trr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trar = ( Chr_r_plus_eps_Ccc.Tra - Chr_r_minus_eps_Ccc.Tra ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trpr = ( Chr_r_plus_eps_Ccc.Trp - Chr_r_minus_eps_Ccc.Trp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tatr = ( Chr_r_plus_eps_Ccc.Tat - Chr_r_minus_eps_Ccc.Tat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tarr = ( Chr_r_plus_eps_Ccc.Tar - Chr_r_minus_eps_Ccc.Tar ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Taar = ( Chr_r_plus_eps_Ccc.Taa - Chr_r_minus_eps_Ccc.Taa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tapr = ( Chr_r_plus_eps_Ccc.Tap - Chr_r_minus_eps_Ccc.Tap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tptr = ( Chr_r_plus_eps_Ccc.Tpt - Chr_r_minus_eps_Ccc.Tpt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tprr = ( Chr_r_plus_eps_Ccc.Tpr - Chr_r_minus_eps_Ccc.Tpr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tpar = ( Chr_r_plus_eps_Ccc.Tpa - Chr_r_minus_eps_Ccc.Tpa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tppr = ( Chr_r_plus_eps_Ccc.Tpp - Chr_r_minus_eps_Ccc.Tpp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rttr = ( Chr_r_plus_eps_Ccc.Rtt - Chr_r_minus_eps_Ccc.Rtt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtrr = ( Chr_r_plus_eps_Ccc.Rtr - Chr_r_minus_eps_Ccc.Rtr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtar = ( Chr_r_plus_eps_Ccc.Rta - Chr_r_minus_eps_Ccc.Rta ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtpr = ( Chr_r_plus_eps_Ccc.Rtp - Chr_r_minus_eps_Ccc.Rtp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrtr = ( Chr_r_plus_eps_Ccc.Rrt - Chr_r_minus_eps_Ccc.Rrt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrrr = ( Chr_r_plus_eps_Ccc.Rrr - Chr_r_minus_eps_Ccc.Rrr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrar = ( Chr_r_plus_eps_Ccc.Rra - Chr_r_minus_eps_Ccc.Rra ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrpr = ( Chr_r_plus_eps_Ccc.Rrp - Chr_r_minus_eps_Ccc.Rrp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ratr = ( Chr_r_plus_eps_Ccc.Rat - Chr_r_minus_eps_Ccc.Rat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rarr = ( Chr_r_plus_eps_Ccc.Rar - Chr_r_minus_eps_Ccc.Rar ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Raar = ( Chr_r_plus_eps_Ccc.Raa - Chr_r_minus_eps_Ccc.Raa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rapr = ( Chr_r_plus_eps_Ccc.Rap - Chr_r_minus_eps_Ccc.Rap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rptr = ( Chr_r_plus_eps_Ccc.Rpt - Chr_r_minus_eps_Ccc.Rpt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rprr = ( Chr_r_plus_eps_Ccc.Rpr - Chr_r_minus_eps_Ccc.Rpr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rpar = ( Chr_r_plus_eps_Ccc.Rpa - Chr_r_minus_eps_Ccc.Rpa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rppr = ( Chr_r_plus_eps_Ccc.Rpp - Chr_r_minus_eps_Ccc.Rpp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Attr = ( Chr_r_plus_eps_Ccc.Att - Chr_r_minus_eps_Ccc.Att ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Atrr = ( Chr_r_plus_eps_Ccc.Atr - Chr_r_minus_eps_Ccc.Atr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Atar = ( Chr_r_plus_eps_Ccc.Ata - Chr_r_minus_eps_Ccc.Ata ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Atpr = ( Chr_r_plus_eps_Ccc.Atp - Chr_r_minus_eps_Ccc.Atp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Artr = ( Chr_r_plus_eps_Ccc.Art - Chr_r_minus_eps_Ccc.Art ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Arrr = ( Chr_r_plus_eps_Ccc.Arr - Chr_r_minus_eps_Ccc.Arr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Arar = ( Chr_r_plus_eps_Ccc.Ara - Chr_r_minus_eps_Ccc.Ara ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Arpr = ( Chr_r_plus_eps_Ccc.Arp - Chr_r_minus_eps_Ccc.Arp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aatr = ( Chr_r_plus_eps_Ccc.Aat - Chr_r_minus_eps_Ccc.Aat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aarr = ( Chr_r_plus_eps_Ccc.Aar - Chr_r_minus_eps_Ccc.Aar ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aaar = ( Chr_r_plus_eps_Ccc.Aaa - Chr_r_minus_eps_Ccc.Aaa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aapr = ( Chr_r_plus_eps_Ccc.Aap - Chr_r_minus_eps_Ccc.Aap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aptr = ( Chr_r_plus_eps_Ccc.Apt - Chr_r_minus_eps_Ccc.Apt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aprr = ( Chr_r_plus_eps_Ccc.Apr - Chr_r_minus_eps_Ccc.Apr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Apar = ( Chr_r_plus_eps_Ccc.Apa - Chr_r_minus_eps_Ccc.Apa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Appr = ( Chr_r_plus_eps_Ccc.App - Chr_r_minus_eps_Ccc.App ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Pttr = ( Chr_r_plus_eps_Ccc.Ptt - Chr_r_minus_eps_Ccc.Ptt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptrr = ( Chr_r_plus_eps_Ccc.Ptr - Chr_r_minus_eps_Ccc.Ptr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptar = ( Chr_r_plus_eps_Ccc.Pta - Chr_r_minus_eps_Ccc.Pta ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptpr = ( Chr_r_plus_eps_Ccc.Ptp - Chr_r_minus_eps_Ccc.Ptp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prtr = ( Chr_r_plus_eps_Ccc.Prt - Chr_r_minus_eps_Ccc.Prt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prrr = ( Chr_r_plus_eps_Ccc.Prr - Chr_r_minus_eps_Ccc.Prr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prar = ( Chr_r_plus_eps_Ccc.Pra - Chr_r_minus_eps_Ccc.Pra ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prpr = ( Chr_r_plus_eps_Ccc.Prp - Chr_r_minus_eps_Ccc.Prp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Patr = ( Chr_r_plus_eps_Ccc.Pat - Chr_r_minus_eps_Ccc.Pat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Parr = ( Chr_r_plus_eps_Ccc.Par - Chr_r_minus_eps_Ccc.Par ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Paar = ( Chr_r_plus_eps_Ccc.Paa - Chr_r_minus_eps_Ccc.Paa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Papr = ( Chr_r_plus_eps_Ccc.Pap - Chr_r_minus_eps_Ccc.Pap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Pptr = ( Chr_r_plus_eps_Ccc.Ppt - Chr_r_minus_eps_Ccc.Ppt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Pprr = ( Chr_r_plus_eps_Ccc.Ppr - Chr_r_minus_eps_Ccc.Ppr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ppar = ( Chr_r_plus_eps_Ccc.Ppa - Chr_r_minus_eps_Ccc.Ppa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Pppr = ( Chr_r_plus_eps_Ccc.Ppp - Chr_r_minus_eps_Ccc.Ppp ) / ( 2 
                    * EPS );

 _d_Chr_Cccc.Ttta = ( Chr_a_plus_eps_Ccc.Ttt - Chr_a_minus_eps_Ccc.Ttt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ttra = ( Chr_a_plus_eps_Ccc.Ttr - Chr_a_minus_eps_Ccc.Ttr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ttaa = ( Chr_a_plus_eps_Ccc.Tta - Chr_a_minus_eps_Ccc.Tta ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ttpa = ( Chr_a_plus_eps_Ccc.Ttp - Chr_a_minus_eps_Ccc.Ttp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trta = ( Chr_a_plus_eps_Ccc.Trt - Chr_a_minus_eps_Ccc.Trt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trra = ( Chr_a_plus_eps_Ccc.Trr - Chr_a_minus_eps_Ccc.Trr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Traa = ( Chr_a_plus_eps_Ccc.Tra - Chr_a_minus_eps_Ccc.Tra ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Trpa = ( Chr_a_plus_eps_Ccc.Trp - Chr_a_minus_eps_Ccc.Trp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tata = ( Chr_a_plus_eps_Ccc.Tat - Chr_a_minus_eps_Ccc.Tat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tara = ( Chr_a_plus_eps_Ccc.Tar - Chr_a_minus_eps_Ccc.Tar ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Taaa = ( Chr_a_plus_eps_Ccc.Taa - Chr_a_minus_eps_Ccc.Taa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tapa = ( Chr_a_plus_eps_Ccc.Tap - Chr_a_minus_eps_Ccc.Tap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tpta = ( Chr_a_plus_eps_Ccc.Tpt - Chr_a_minus_eps_Ccc.Tpt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tpra = ( Chr_a_plus_eps_Ccc.Tpr - Chr_a_minus_eps_Ccc.Tpr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tpaa = ( Chr_a_plus_eps_Ccc.Tpa - Chr_a_minus_eps_Ccc.Tpa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Tppa = ( Chr_a_plus_eps_Ccc.Tpp - Chr_a_minus_eps_Ccc.Tpp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtta = ( Chr_a_plus_eps_Ccc.Rtt - Chr_a_minus_eps_Ccc.Rtt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtra = ( Chr_a_plus_eps_Ccc.Rtr - Chr_a_minus_eps_Ccc.Rtr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtaa = ( Chr_a_plus_eps_Ccc.Rta - Chr_a_minus_eps_Ccc.Rta ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rtpa = ( Chr_a_plus_eps_Ccc.Rtp - Chr_a_minus_eps_Ccc.Rtp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrta = ( Chr_a_plus_eps_Ccc.Rrt - Chr_a_minus_eps_Ccc.Rrt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrra = ( Chr_a_plus_eps_Ccc.Rrr - Chr_a_minus_eps_Ccc.Rrr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rraa = ( Chr_a_plus_eps_Ccc.Rra - Chr_a_minus_eps_Ccc.Rra ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rrpa = ( Chr_a_plus_eps_Ccc.Rrp - Chr_a_minus_eps_Ccc.Rrp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rata = ( Chr_a_plus_eps_Ccc.Rat - Chr_a_minus_eps_Ccc.Rat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rara = ( Chr_a_plus_eps_Ccc.Rar - Chr_a_minus_eps_Ccc.Rar ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Raaa = ( Chr_a_plus_eps_Ccc.Raa - Chr_a_minus_eps_Ccc.Raa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rapa = ( Chr_a_plus_eps_Ccc.Rap - Chr_a_minus_eps_Ccc.Rap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rpta = ( Chr_a_plus_eps_Ccc.Rpt - Chr_a_minus_eps_Ccc.Rpt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rpra = ( Chr_a_plus_eps_Ccc.Rpr - Chr_a_minus_eps_Ccc.Rpr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rpaa = ( Chr_a_plus_eps_Ccc.Rpa - Chr_a_minus_eps_Ccc.Rpa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Rppa = ( Chr_a_plus_eps_Ccc.Rpp - Chr_a_minus_eps_Ccc.Rpp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Atta = ( Chr_a_plus_eps_Ccc.Att - Chr_a_minus_eps_Ccc.Att ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Atra = ( Chr_a_plus_eps_Ccc.Atr - Chr_a_minus_eps_Ccc.Atr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ataa = ( Chr_a_plus_eps_Ccc.Ata - Chr_a_minus_eps_Ccc.Ata ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Atpa = ( Chr_a_plus_eps_Ccc.Atp - Chr_a_minus_eps_Ccc.Atp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Arta = ( Chr_a_plus_eps_Ccc.Art - Chr_a_minus_eps_Ccc.Art ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Arra = ( Chr_a_plus_eps_Ccc.Arr - Chr_a_minus_eps_Ccc.Arr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Araa = ( Chr_a_plus_eps_Ccc.Ara - Chr_a_minus_eps_Ccc.Ara ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Arpa = ( Chr_a_plus_eps_Ccc.Arp - Chr_a_minus_eps_Ccc.Arp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aata = ( Chr_a_plus_eps_Ccc.Aat - Chr_a_minus_eps_Ccc.Aat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aara = ( Chr_a_plus_eps_Ccc.Aar - Chr_a_minus_eps_Ccc.Aar ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aaaa = ( Chr_a_plus_eps_Ccc.Aaa - Chr_a_minus_eps_Ccc.Aaa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Aapa = ( Chr_a_plus_eps_Ccc.Aap - Chr_a_minus_eps_Ccc.Aap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Apta = ( Chr_a_plus_eps_Ccc.Apt - Chr_a_minus_eps_Ccc.Apt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Apra = ( Chr_a_plus_eps_Ccc.Apr - Chr_a_minus_eps_Ccc.Apr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Apaa = ( Chr_a_plus_eps_Ccc.Apa - Chr_a_minus_eps_Ccc.Apa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Appa = ( Chr_a_plus_eps_Ccc.App - Chr_a_minus_eps_Ccc.App ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptta = ( Chr_a_plus_eps_Ccc.Ptt - Chr_a_minus_eps_Ccc.Ptt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptra = ( Chr_a_plus_eps_Ccc.Ptr - Chr_a_minus_eps_Ccc.Ptr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptaa = ( Chr_a_plus_eps_Ccc.Pta - Chr_a_minus_eps_Ccc.Pta ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ptpa = ( Chr_a_plus_eps_Ccc.Ptp - Chr_a_minus_eps_Ccc.Ptp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prta = ( Chr_a_plus_eps_Ccc.Prt - Chr_a_minus_eps_Ccc.Prt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prra = ( Chr_a_plus_eps_Ccc.Prr - Chr_a_minus_eps_Ccc.Prr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Praa = ( Chr_a_plus_eps_Ccc.Pra - Chr_a_minus_eps_Ccc.Pra ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Prpa = ( Chr_a_plus_eps_Ccc.Prp - Chr_a_minus_eps_Ccc.Prp ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Pata = ( Chr_a_plus_eps_Ccc.Pat - Chr_a_minus_eps_Ccc.Pat ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Para = ( Chr_a_plus_eps_Ccc.Par - Chr_a_minus_eps_Ccc.Par ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Paaa = ( Chr_a_plus_eps_Ccc.Paa - Chr_a_minus_eps_Ccc.Paa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Papa = ( Chr_a_plus_eps_Ccc.Pap - Chr_a_minus_eps_Ccc.Pap ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ppta = ( Chr_a_plus_eps_Ccc.Ppt - Chr_a_minus_eps_Ccc.Ppt ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ppra = ( Chr_a_plus_eps_Ccc.Ppr - Chr_a_minus_eps_Ccc.Ppr ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Ppaa = ( Chr_a_plus_eps_Ccc.Ppa - Chr_a_minus_eps_Ccc.Ppa ) / ( 2 
                    * EPS );
 _d_Chr_Cccc.Pppa = ( Chr_a_plus_eps_Ccc.Ppp - Chr_a_minus_eps_Ccc.Ppp ) / ( 2 
                    * EPS );

 _d_Chr_Cccc.Tttp = 0.0;
 _d_Chr_Cccc.Ttrp = 0.0;
 _d_Chr_Cccc.Ttap = 0.0;
 _d_Chr_Cccc.Ttpp = 0.0;
 _d_Chr_Cccc.Trtp = 0.0;
 _d_Chr_Cccc.Trrp = 0.0;
 _d_Chr_Cccc.Trap = 0.0;
 _d_Chr_Cccc.Trpp = 0.0;
 _d_Chr_Cccc.Tatp = 0.0;
 _d_Chr_Cccc.Tarp = 0.0;
 _d_Chr_Cccc.Taap = 0.0;
 _d_Chr_Cccc.Tapp = 0.0;
 _d_Chr_Cccc.Tptp = 0.0;
 _d_Chr_Cccc.Tprp = 0.0;
 _d_Chr_Cccc.Tpap = 0.0;
 _d_Chr_Cccc.Tppp = 0.0;
 _d_Chr_Cccc.Rttp = 0.0;
 _d_Chr_Cccc.Rtrp = 0.0;
 _d_Chr_Cccc.Rtap = 0.0;
 _d_Chr_Cccc.Rtpp = 0.0;
 _d_Chr_Cccc.Rrtp = 0.0;
 _d_Chr_Cccc.Rrrp = 0.0;
 _d_Chr_Cccc.Rrap = 0.0;
 _d_Chr_Cccc.Rrpp = 0.0;
 _d_Chr_Cccc.Ratp = 0.0;
 _d_Chr_Cccc.Rarp = 0.0;
 _d_Chr_Cccc.Raap = 0.0;
 _d_Chr_Cccc.Rapp = 0.0;
 _d_Chr_Cccc.Rptp = 0.0;
 _d_Chr_Cccc.Rprp = 0.0;
 _d_Chr_Cccc.Rpap = 0.0;
 _d_Chr_Cccc.Rppp = 0.0;
 _d_Chr_Cccc.Attp = 0.0;
 _d_Chr_Cccc.Atrp = 0.0;
 _d_Chr_Cccc.Atap = 0.0;
 _d_Chr_Cccc.Atpp = 0.0;
 _d_Chr_Cccc.Artp = 0.0;
 _d_Chr_Cccc.Arrp = 0.0;
 _d_Chr_Cccc.Arap = 0.0;
 _d_Chr_Cccc.Arpp = 0.0;
 _d_Chr_Cccc.Aatp = 0.0;
 _d_Chr_Cccc.Aarp = 0.0;
 _d_Chr_Cccc.Aaap = 0.0;
 _d_Chr_Cccc.Aapp = 0.0;
 _d_Chr_Cccc.Aptp = 0.0;
 _d_Chr_Cccc.Aprp = 0.0;
 _d_Chr_Cccc.Apap = 0.0;
 _d_Chr_Cccc.Appp = 0.0;
 _d_Chr_Cccc.Pttp = 0.0;
 _d_Chr_Cccc.Ptrp = 0.0;
 _d_Chr_Cccc.Ptap = 0.0;
 _d_Chr_Cccc.Ptpp = 0.0;
 _d_Chr_Cccc.Prtp = 0.0;
 _d_Chr_Cccc.Prrp = 0.0;
 _d_Chr_Cccc.Prap = 0.0;
 _d_Chr_Cccc.Prpp = 0.0;
 _d_Chr_Cccc.Patp = 0.0;
 _d_Chr_Cccc.Parp = 0.0;
 _d_Chr_Cccc.Paap = 0.0;
 _d_Chr_Cccc.Papp = 0.0;
 _d_Chr_Cccc.Pptp = 0.0;
 _d_Chr_Cccc.Pprp = 0.0;
 _d_Chr_Cccc.Ppap = 0.0;
 _d_Chr_Cccc.Pppp = 0.0;

 R_Cccc.Tttt = _d_Chr_Cccc.Tttt - _d_Chr_Cccc.Tttt + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Ttt + Chr_Ccc.Trt * Chr_Ccc.Rtt + Chr_Ccc.Tat 
               * Chr_Ccc.Att + Chr_Ccc.Tpt * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Ttt + Chr_Ccc.Trt * Chr_Ccc.Rtt + Chr_Ccc.Tat 
               * Chr_Ccc.Att + Chr_Ccc.Tpt * Chr_Ccc.Ptt );
 R_Cccc.Tttr = _d_Chr_Cccc.Ttrt - _d_Chr_Cccc.Tttr + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Ttr + Chr_Ccc.Trt * Chr_Ccc.Rtr + Chr_Ccc.Tat 
               * Chr_Ccc.Atr + Chr_Ccc.Tpt * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Ttt + Chr_Ccc.Trr * Chr_Ccc.Rtt + Chr_Ccc.Tar 
               * Chr_Ccc.Att + Chr_Ccc.Tpr * Chr_Ccc.Ptt );
 R_Cccc.Ttta = _d_Chr_Cccc.Ttat - _d_Chr_Cccc.Ttta + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tta + Chr_Ccc.Trt * Chr_Ccc.Rta + Chr_Ccc.Tat 
               * Chr_Ccc.Ata + Chr_Ccc.Tpt * Chr_Ccc.Pta ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Ttt + Chr_Ccc.Tra * Chr_Ccc.Rtt + Chr_Ccc.Taa 
               * Chr_Ccc.Att + Chr_Ccc.Tpa * Chr_Ccc.Ptt );
 R_Cccc.Tttp = _d_Chr_Cccc.Ttpt - _d_Chr_Cccc.Tttp + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Ttp + Chr_Ccc.Trt * Chr_Ccc.Rtp + Chr_Ccc.Tat 
               * Chr_Ccc.Atp + Chr_Ccc.Tpt * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Ttt + Chr_Ccc.Trp * Chr_Ccc.Rtt + Chr_Ccc.Tap 
               * Chr_Ccc.Att + Chr_Ccc.Tpp * Chr_Ccc.Ptt );
 R_Cccc.Ttrt = _d_Chr_Cccc.Tttr - _d_Chr_Cccc.Ttrt + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Ttt + Chr_Ccc.Trr * Chr_Ccc.Rtt + Chr_Ccc.Tar 
               * Chr_Ccc.Att + Chr_Ccc.Tpr * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Ttr + Chr_Ccc.Trt * Chr_Ccc.Rtr + Chr_Ccc.Tat 
               * Chr_Ccc.Atr + Chr_Ccc.Tpt * Chr_Ccc.Ptr );
 R_Cccc.Ttrr = _d_Chr_Cccc.Ttrr - _d_Chr_Cccc.Ttrr + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Ttr + Chr_Ccc.Trr * Chr_Ccc.Rtr + Chr_Ccc.Tar 
               * Chr_Ccc.Atr + Chr_Ccc.Tpr * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Ttr + Chr_Ccc.Trr * Chr_Ccc.Rtr + Chr_Ccc.Tar 
               * Chr_Ccc.Atr + Chr_Ccc.Tpr * Chr_Ccc.Ptr );
 R_Cccc.Ttra = _d_Chr_Cccc.Ttar - _d_Chr_Cccc.Ttra + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tta + Chr_Ccc.Trr * Chr_Ccc.Rta + Chr_Ccc.Tar 
               * Chr_Ccc.Ata + Chr_Ccc.Tpr * Chr_Ccc.Pta ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Ttr + Chr_Ccc.Tra * Chr_Ccc.Rtr + Chr_Ccc.Taa 
               * Chr_Ccc.Atr + Chr_Ccc.Tpa * Chr_Ccc.Ptr );
 R_Cccc.Ttrp = _d_Chr_Cccc.Ttpr - _d_Chr_Cccc.Ttrp + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Ttp + Chr_Ccc.Trr * Chr_Ccc.Rtp + Chr_Ccc.Tar 
               * Chr_Ccc.Atp + Chr_Ccc.Tpr * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Ttr + Chr_Ccc.Trp * Chr_Ccc.Rtr + Chr_Ccc.Tap 
               * Chr_Ccc.Atr + Chr_Ccc.Tpp * Chr_Ccc.Ptr );
 R_Cccc.Ttat = _d_Chr_Cccc.Ttta - _d_Chr_Cccc.Ttat + ( Chr_Ccc.Tta 
               * Chr_Ccc.Ttt + Chr_Ccc.Tra * Chr_Ccc.Rtt + Chr_Ccc.Taa 
               * Chr_Ccc.Att + Chr_Ccc.Tpa * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tta + Chr_Ccc.Trt * Chr_Ccc.Rta + Chr_Ccc.Tat 
               * Chr_Ccc.Ata + Chr_Ccc.Tpt * Chr_Ccc.Pta );
 R_Cccc.Ttar = _d_Chr_Cccc.Ttra - _d_Chr_Cccc.Ttar + ( Chr_Ccc.Tta 
               * Chr_Ccc.Ttr + Chr_Ccc.Tra * Chr_Ccc.Rtr + Chr_Ccc.Taa 
               * Chr_Ccc.Atr + Chr_Ccc.Tpa * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tta + Chr_Ccc.Trr * Chr_Ccc.Rta + Chr_Ccc.Tar 
               * Chr_Ccc.Ata + Chr_Ccc.Tpr * Chr_Ccc.Pta );
 R_Cccc.Ttaa = _d_Chr_Cccc.Ttaa - _d_Chr_Cccc.Ttaa + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tta + Chr_Ccc.Tra * Chr_Ccc.Rta + Chr_Ccc.Taa 
               * Chr_Ccc.Ata + Chr_Ccc.Tpa * Chr_Ccc.Pta ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tta + Chr_Ccc.Tra * Chr_Ccc.Rta + Chr_Ccc.Taa 
               * Chr_Ccc.Ata + Chr_Ccc.Tpa * Chr_Ccc.Pta );
 R_Cccc.Ttap = _d_Chr_Cccc.Ttpa - _d_Chr_Cccc.Ttap + ( Chr_Ccc.Tta 
               * Chr_Ccc.Ttp + Chr_Ccc.Tra * Chr_Ccc.Rtp + Chr_Ccc.Taa 
               * Chr_Ccc.Atp + Chr_Ccc.Tpa * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tta + Chr_Ccc.Trp * Chr_Ccc.Rta + Chr_Ccc.Tap 
               * Chr_Ccc.Ata + Chr_Ccc.Tpp * Chr_Ccc.Pta );
 R_Cccc.Ttpt = _d_Chr_Cccc.Tttp - _d_Chr_Cccc.Ttpt + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Ttt + Chr_Ccc.Trp * Chr_Ccc.Rtt + Chr_Ccc.Tap 
               * Chr_Ccc.Att + Chr_Ccc.Tpp * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Ttp + Chr_Ccc.Trt * Chr_Ccc.Rtp + Chr_Ccc.Tat 
               * Chr_Ccc.Atp + Chr_Ccc.Tpt * Chr_Ccc.Ptp );
 R_Cccc.Ttpr = _d_Chr_Cccc.Ttrp - _d_Chr_Cccc.Ttpr + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Ttr + Chr_Ccc.Trp * Chr_Ccc.Rtr + Chr_Ccc.Tap 
               * Chr_Ccc.Atr + Chr_Ccc.Tpp * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Ttp + Chr_Ccc.Trr * Chr_Ccc.Rtp + Chr_Ccc.Tar 
               * Chr_Ccc.Atp + Chr_Ccc.Tpr * Chr_Ccc.Ptp );
 R_Cccc.Ttpa = _d_Chr_Cccc.Ttap - _d_Chr_Cccc.Ttpa + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tta + Chr_Ccc.Trp * Chr_Ccc.Rta + Chr_Ccc.Tap 
               * Chr_Ccc.Ata + Chr_Ccc.Tpp * Chr_Ccc.Pta ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Ttp + Chr_Ccc.Tra * Chr_Ccc.Rtp + Chr_Ccc.Taa 
               * Chr_Ccc.Atp + Chr_Ccc.Tpa * Chr_Ccc.Ptp );
 R_Cccc.Ttpp = _d_Chr_Cccc.Ttpp - _d_Chr_Cccc.Ttpp + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Ttp + Chr_Ccc.Trp * Chr_Ccc.Rtp + Chr_Ccc.Tap 
               * Chr_Ccc.Atp + Chr_Ccc.Tpp * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Ttp + Chr_Ccc.Trp * Chr_Ccc.Rtp + Chr_Ccc.Tap 
               * Chr_Ccc.Atp + Chr_Ccc.Tpp * Chr_Ccc.Ptp );
 R_Cccc.Trtt = _d_Chr_Cccc.Trtt - _d_Chr_Cccc.Trtt + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Trt + Chr_Ccc.Trt * Chr_Ccc.Rrt + Chr_Ccc.Tat 
               * Chr_Ccc.Art + Chr_Ccc.Tpt * Chr_Ccc.Prt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Trt + Chr_Ccc.Trt * Chr_Ccc.Rrt + Chr_Ccc.Tat 
               * Chr_Ccc.Art + Chr_Ccc.Tpt * Chr_Ccc.Prt );
 R_Cccc.Trtr = _d_Chr_Cccc.Trrt - _d_Chr_Cccc.Trtr + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Trr + Chr_Ccc.Trt * Chr_Ccc.Rrr + Chr_Ccc.Tat 
               * Chr_Ccc.Arr + Chr_Ccc.Tpt * Chr_Ccc.Prr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Trt + Chr_Ccc.Trr * Chr_Ccc.Rrt + Chr_Ccc.Tar 
               * Chr_Ccc.Art + Chr_Ccc.Tpr * Chr_Ccc.Prt );
 R_Cccc.Trta = _d_Chr_Cccc.Trat - _d_Chr_Cccc.Trta + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tra + Chr_Ccc.Trt * Chr_Ccc.Rra + Chr_Ccc.Tat 
               * Chr_Ccc.Ara + Chr_Ccc.Tpt * Chr_Ccc.Pra ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Trt + Chr_Ccc.Tra * Chr_Ccc.Rrt + Chr_Ccc.Taa 
               * Chr_Ccc.Art + Chr_Ccc.Tpa * Chr_Ccc.Prt );
 R_Cccc.Trtp = _d_Chr_Cccc.Trpt - _d_Chr_Cccc.Trtp + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Trp + Chr_Ccc.Trt * Chr_Ccc.Rrp + Chr_Ccc.Tat 
               * Chr_Ccc.Arp + Chr_Ccc.Tpt * Chr_Ccc.Prp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Trt + Chr_Ccc.Trp * Chr_Ccc.Rrt + Chr_Ccc.Tap 
               * Chr_Ccc.Art + Chr_Ccc.Tpp * Chr_Ccc.Prt );
 R_Cccc.Trrt = _d_Chr_Cccc.Trtr - _d_Chr_Cccc.Trrt + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Trt + Chr_Ccc.Trr * Chr_Ccc.Rrt + Chr_Ccc.Tar 
               * Chr_Ccc.Art + Chr_Ccc.Tpr * Chr_Ccc.Prt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Trr + Chr_Ccc.Trt * Chr_Ccc.Rrr + Chr_Ccc.Tat 
               * Chr_Ccc.Arr + Chr_Ccc.Tpt * Chr_Ccc.Prr );
 R_Cccc.Trrr = _d_Chr_Cccc.Trrr - _d_Chr_Cccc.Trrr + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Trr + Chr_Ccc.Trr * Chr_Ccc.Rrr + Chr_Ccc.Tar 
               * Chr_Ccc.Arr + Chr_Ccc.Tpr * Chr_Ccc.Prr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Trr + Chr_Ccc.Trr * Chr_Ccc.Rrr + Chr_Ccc.Tar 
               * Chr_Ccc.Arr + Chr_Ccc.Tpr * Chr_Ccc.Prr );
 R_Cccc.Trra = _d_Chr_Cccc.Trar - _d_Chr_Cccc.Trra + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tra + Chr_Ccc.Trr * Chr_Ccc.Rra + Chr_Ccc.Tar 
               * Chr_Ccc.Ara + Chr_Ccc.Tpr * Chr_Ccc.Pra ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Trr + Chr_Ccc.Tra * Chr_Ccc.Rrr + Chr_Ccc.Taa 
               * Chr_Ccc.Arr + Chr_Ccc.Tpa * Chr_Ccc.Prr );
 R_Cccc.Trrp = _d_Chr_Cccc.Trpr - _d_Chr_Cccc.Trrp + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Trp + Chr_Ccc.Trr * Chr_Ccc.Rrp + Chr_Ccc.Tar 
               * Chr_Ccc.Arp + Chr_Ccc.Tpr * Chr_Ccc.Prp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Trr + Chr_Ccc.Trp * Chr_Ccc.Rrr + Chr_Ccc.Tap 
               * Chr_Ccc.Arr + Chr_Ccc.Tpp * Chr_Ccc.Prr );
 R_Cccc.Trat = _d_Chr_Cccc.Trta - _d_Chr_Cccc.Trat + ( Chr_Ccc.Tta 
               * Chr_Ccc.Trt + Chr_Ccc.Tra * Chr_Ccc.Rrt + Chr_Ccc.Taa 
               * Chr_Ccc.Art + Chr_Ccc.Tpa * Chr_Ccc.Prt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tra + Chr_Ccc.Trt * Chr_Ccc.Rra + Chr_Ccc.Tat 
               * Chr_Ccc.Ara + Chr_Ccc.Tpt * Chr_Ccc.Pra );
 R_Cccc.Trar = _d_Chr_Cccc.Trra - _d_Chr_Cccc.Trar + ( Chr_Ccc.Tta 
               * Chr_Ccc.Trr + Chr_Ccc.Tra * Chr_Ccc.Rrr + Chr_Ccc.Taa 
               * Chr_Ccc.Arr + Chr_Ccc.Tpa * Chr_Ccc.Prr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tra + Chr_Ccc.Trr * Chr_Ccc.Rra + Chr_Ccc.Tar 
               * Chr_Ccc.Ara + Chr_Ccc.Tpr * Chr_Ccc.Pra );
 R_Cccc.Traa = _d_Chr_Cccc.Traa - _d_Chr_Cccc.Traa + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tra + Chr_Ccc.Tra * Chr_Ccc.Rra + Chr_Ccc.Taa 
               * Chr_Ccc.Ara + Chr_Ccc.Tpa * Chr_Ccc.Pra ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tra + Chr_Ccc.Tra * Chr_Ccc.Rra + Chr_Ccc.Taa 
               * Chr_Ccc.Ara + Chr_Ccc.Tpa * Chr_Ccc.Pra );
 R_Cccc.Trap = _d_Chr_Cccc.Trpa - _d_Chr_Cccc.Trap + ( Chr_Ccc.Tta 
               * Chr_Ccc.Trp + Chr_Ccc.Tra * Chr_Ccc.Rrp + Chr_Ccc.Taa 
               * Chr_Ccc.Arp + Chr_Ccc.Tpa * Chr_Ccc.Prp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tra + Chr_Ccc.Trp * Chr_Ccc.Rra + Chr_Ccc.Tap 
               * Chr_Ccc.Ara + Chr_Ccc.Tpp * Chr_Ccc.Pra );
 R_Cccc.Trpt = _d_Chr_Cccc.Trtp - _d_Chr_Cccc.Trpt + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Trt + Chr_Ccc.Trp * Chr_Ccc.Rrt + Chr_Ccc.Tap 
               * Chr_Ccc.Art + Chr_Ccc.Tpp * Chr_Ccc.Prt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Trp + Chr_Ccc.Trt * Chr_Ccc.Rrp + Chr_Ccc.Tat 
               * Chr_Ccc.Arp + Chr_Ccc.Tpt * Chr_Ccc.Prp );
 R_Cccc.Trpr = _d_Chr_Cccc.Trrp - _d_Chr_Cccc.Trpr + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Trr + Chr_Ccc.Trp * Chr_Ccc.Rrr + Chr_Ccc.Tap 
               * Chr_Ccc.Arr + Chr_Ccc.Tpp * Chr_Ccc.Prr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Trp + Chr_Ccc.Trr * Chr_Ccc.Rrp + Chr_Ccc.Tar 
               * Chr_Ccc.Arp + Chr_Ccc.Tpr * Chr_Ccc.Prp );
 R_Cccc.Trpa = _d_Chr_Cccc.Trap - _d_Chr_Cccc.Trpa + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tra + Chr_Ccc.Trp * Chr_Ccc.Rra + Chr_Ccc.Tap 
               * Chr_Ccc.Ara + Chr_Ccc.Tpp * Chr_Ccc.Pra ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Trp + Chr_Ccc.Tra * Chr_Ccc.Rrp + Chr_Ccc.Taa 
               * Chr_Ccc.Arp + Chr_Ccc.Tpa * Chr_Ccc.Prp );
 R_Cccc.Trpp = _d_Chr_Cccc.Trpp - _d_Chr_Cccc.Trpp + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Trp + Chr_Ccc.Trp * Chr_Ccc.Rrp + Chr_Ccc.Tap 
               * Chr_Ccc.Arp + Chr_Ccc.Tpp * Chr_Ccc.Prp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Trp + Chr_Ccc.Trp * Chr_Ccc.Rrp + Chr_Ccc.Tap 
               * Chr_Ccc.Arp + Chr_Ccc.Tpp * Chr_Ccc.Prp );
 R_Cccc.Tatt = _d_Chr_Cccc.Tatt - _d_Chr_Cccc.Tatt + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tat + Chr_Ccc.Trt * Chr_Ccc.Rat + Chr_Ccc.Tat 
               * Chr_Ccc.Aat + Chr_Ccc.Tpt * Chr_Ccc.Pat ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tat + Chr_Ccc.Trt * Chr_Ccc.Rat + Chr_Ccc.Tat 
               * Chr_Ccc.Aat + Chr_Ccc.Tpt * Chr_Ccc.Pat );
 R_Cccc.Tatr = _d_Chr_Cccc.Tart - _d_Chr_Cccc.Tatr + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tar + Chr_Ccc.Trt * Chr_Ccc.Rar + Chr_Ccc.Tat 
               * Chr_Ccc.Aar + Chr_Ccc.Tpt * Chr_Ccc.Par ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tat + Chr_Ccc.Trr * Chr_Ccc.Rat + Chr_Ccc.Tar 
               * Chr_Ccc.Aat + Chr_Ccc.Tpr * Chr_Ccc.Pat );
 R_Cccc.Tata = _d_Chr_Cccc.Taat - _d_Chr_Cccc.Tata + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Taa + Chr_Ccc.Trt * Chr_Ccc.Raa + Chr_Ccc.Tat 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpt * Chr_Ccc.Paa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tat + Chr_Ccc.Tra * Chr_Ccc.Rat + Chr_Ccc.Taa 
               * Chr_Ccc.Aat + Chr_Ccc.Tpa * Chr_Ccc.Pat );
 R_Cccc.Tatp = _d_Chr_Cccc.Tapt - _d_Chr_Cccc.Tatp + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tap + Chr_Ccc.Trt * Chr_Ccc.Rap + Chr_Ccc.Tat 
               * Chr_Ccc.Aap + Chr_Ccc.Tpt * Chr_Ccc.Pap ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tat + Chr_Ccc.Trp * Chr_Ccc.Rat + Chr_Ccc.Tap 
               * Chr_Ccc.Aat + Chr_Ccc.Tpp * Chr_Ccc.Pat );
 R_Cccc.Tart = _d_Chr_Cccc.Tatr - _d_Chr_Cccc.Tart + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tat + Chr_Ccc.Trr * Chr_Ccc.Rat + Chr_Ccc.Tar 
               * Chr_Ccc.Aat + Chr_Ccc.Tpr * Chr_Ccc.Pat ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tar + Chr_Ccc.Trt * Chr_Ccc.Rar + Chr_Ccc.Tat 
               * Chr_Ccc.Aar + Chr_Ccc.Tpt * Chr_Ccc.Par );
 R_Cccc.Tarr = _d_Chr_Cccc.Tarr - _d_Chr_Cccc.Tarr + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tar + Chr_Ccc.Trr * Chr_Ccc.Rar + Chr_Ccc.Tar 
               * Chr_Ccc.Aar + Chr_Ccc.Tpr * Chr_Ccc.Par ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tar + Chr_Ccc.Trr * Chr_Ccc.Rar + Chr_Ccc.Tar 
               * Chr_Ccc.Aar + Chr_Ccc.Tpr * Chr_Ccc.Par );
 R_Cccc.Tara = _d_Chr_Cccc.Taar - _d_Chr_Cccc.Tara + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Taa + Chr_Ccc.Trr * Chr_Ccc.Raa + Chr_Ccc.Tar 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpr * Chr_Ccc.Paa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tar + Chr_Ccc.Tra * Chr_Ccc.Rar + Chr_Ccc.Taa 
               * Chr_Ccc.Aar + Chr_Ccc.Tpa * Chr_Ccc.Par );
 R_Cccc.Tarp = _d_Chr_Cccc.Tapr - _d_Chr_Cccc.Tarp + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tap + Chr_Ccc.Trr * Chr_Ccc.Rap + Chr_Ccc.Tar 
               * Chr_Ccc.Aap + Chr_Ccc.Tpr * Chr_Ccc.Pap ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tar + Chr_Ccc.Trp * Chr_Ccc.Rar + Chr_Ccc.Tap 
               * Chr_Ccc.Aar + Chr_Ccc.Tpp * Chr_Ccc.Par );
 R_Cccc.Taat = _d_Chr_Cccc.Tata - _d_Chr_Cccc.Taat + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tat + Chr_Ccc.Tra * Chr_Ccc.Rat + Chr_Ccc.Taa 
               * Chr_Ccc.Aat + Chr_Ccc.Tpa * Chr_Ccc.Pat ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Taa + Chr_Ccc.Trt * Chr_Ccc.Raa + Chr_Ccc.Tat 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpt * Chr_Ccc.Paa );
 R_Cccc.Taar = _d_Chr_Cccc.Tara - _d_Chr_Cccc.Taar + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tar + Chr_Ccc.Tra * Chr_Ccc.Rar + Chr_Ccc.Taa 
               * Chr_Ccc.Aar + Chr_Ccc.Tpa * Chr_Ccc.Par ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Taa + Chr_Ccc.Trr * Chr_Ccc.Raa + Chr_Ccc.Tar 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpr * Chr_Ccc.Paa );
 R_Cccc.Taaa = _d_Chr_Cccc.Taaa - _d_Chr_Cccc.Taaa + ( Chr_Ccc.Tta 
               * Chr_Ccc.Taa + Chr_Ccc.Tra * Chr_Ccc.Raa + Chr_Ccc.Taa 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpa * Chr_Ccc.Paa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Taa + Chr_Ccc.Tra * Chr_Ccc.Raa + Chr_Ccc.Taa 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpa * Chr_Ccc.Paa );
 R_Cccc.Taap = _d_Chr_Cccc.Tapa - _d_Chr_Cccc.Taap + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tap + Chr_Ccc.Tra * Chr_Ccc.Rap + Chr_Ccc.Taa 
               * Chr_Ccc.Aap + Chr_Ccc.Tpa * Chr_Ccc.Pap ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Taa + Chr_Ccc.Trp * Chr_Ccc.Raa + Chr_Ccc.Tap 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpp * Chr_Ccc.Paa );
 R_Cccc.Tapt = _d_Chr_Cccc.Tatp - _d_Chr_Cccc.Tapt + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tat + Chr_Ccc.Trp * Chr_Ccc.Rat + Chr_Ccc.Tap 
               * Chr_Ccc.Aat + Chr_Ccc.Tpp * Chr_Ccc.Pat ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tap + Chr_Ccc.Trt * Chr_Ccc.Rap + Chr_Ccc.Tat 
               * Chr_Ccc.Aap + Chr_Ccc.Tpt * Chr_Ccc.Pap );
 R_Cccc.Tapr = _d_Chr_Cccc.Tarp - _d_Chr_Cccc.Tapr + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tar + Chr_Ccc.Trp * Chr_Ccc.Rar + Chr_Ccc.Tap 
               * Chr_Ccc.Aar + Chr_Ccc.Tpp * Chr_Ccc.Par ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tap + Chr_Ccc.Trr * Chr_Ccc.Rap + Chr_Ccc.Tar 
               * Chr_Ccc.Aap + Chr_Ccc.Tpr * Chr_Ccc.Pap );
 R_Cccc.Tapa = _d_Chr_Cccc.Taap - _d_Chr_Cccc.Tapa + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Taa + Chr_Ccc.Trp * Chr_Ccc.Raa + Chr_Ccc.Tap 
               * Chr_Ccc.Aaa + Chr_Ccc.Tpp * Chr_Ccc.Paa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tap + Chr_Ccc.Tra * Chr_Ccc.Rap + Chr_Ccc.Taa 
               * Chr_Ccc.Aap + Chr_Ccc.Tpa * Chr_Ccc.Pap );
 R_Cccc.Tapp = _d_Chr_Cccc.Tapp - _d_Chr_Cccc.Tapp + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tap + Chr_Ccc.Trp * Chr_Ccc.Rap + Chr_Ccc.Tap 
               * Chr_Ccc.Aap + Chr_Ccc.Tpp * Chr_Ccc.Pap ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tap + Chr_Ccc.Trp * Chr_Ccc.Rap + Chr_Ccc.Tap 
               * Chr_Ccc.Aap + Chr_Ccc.Tpp * Chr_Ccc.Pap );
 R_Cccc.Tptt = _d_Chr_Cccc.Tptt - _d_Chr_Cccc.Tptt + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpt + Chr_Ccc.Trt * Chr_Ccc.Rpt + Chr_Ccc.Tat 
               * Chr_Ccc.Apt + Chr_Ccc.Tpt * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpt + Chr_Ccc.Trt * Chr_Ccc.Rpt + Chr_Ccc.Tat 
               * Chr_Ccc.Apt + Chr_Ccc.Tpt * Chr_Ccc.Ppt );
 R_Cccc.Tptr = _d_Chr_Cccc.Tprt - _d_Chr_Cccc.Tptr + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpr + Chr_Ccc.Trt * Chr_Ccc.Rpr + Chr_Ccc.Tat 
               * Chr_Ccc.Apr + Chr_Ccc.Tpt * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpt + Chr_Ccc.Trr * Chr_Ccc.Rpt + Chr_Ccc.Tar 
               * Chr_Ccc.Apt + Chr_Ccc.Tpr * Chr_Ccc.Ppt );
 R_Cccc.Tpta = _d_Chr_Cccc.Tpat - _d_Chr_Cccc.Tpta + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpa + Chr_Ccc.Trt * Chr_Ccc.Rpa + Chr_Ccc.Tat 
               * Chr_Ccc.Apa + Chr_Ccc.Tpt * Chr_Ccc.Ppa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpt + Chr_Ccc.Tra * Chr_Ccc.Rpt + Chr_Ccc.Taa 
               * Chr_Ccc.Apt + Chr_Ccc.Tpa * Chr_Ccc.Ppt );
 R_Cccc.Tptp = _d_Chr_Cccc.Tppt - _d_Chr_Cccc.Tptp + ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpp + Chr_Ccc.Trt * Chr_Ccc.Rpp + Chr_Ccc.Tat 
               * Chr_Ccc.App + Chr_Ccc.Tpt * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpt + Chr_Ccc.Trp * Chr_Ccc.Rpt + Chr_Ccc.Tap 
               * Chr_Ccc.Apt + Chr_Ccc.Tpp * Chr_Ccc.Ppt );
 R_Cccc.Tprt = _d_Chr_Cccc.Tptr - _d_Chr_Cccc.Tprt + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpt + Chr_Ccc.Trr * Chr_Ccc.Rpt + Chr_Ccc.Tar 
               * Chr_Ccc.Apt + Chr_Ccc.Tpr * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpr + Chr_Ccc.Trt * Chr_Ccc.Rpr + Chr_Ccc.Tat 
               * Chr_Ccc.Apr + Chr_Ccc.Tpt * Chr_Ccc.Ppr );
 R_Cccc.Tprr = _d_Chr_Cccc.Tprr - _d_Chr_Cccc.Tprr + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpr + Chr_Ccc.Trr * Chr_Ccc.Rpr + Chr_Ccc.Tar 
               * Chr_Ccc.Apr + Chr_Ccc.Tpr * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpr + Chr_Ccc.Trr * Chr_Ccc.Rpr + Chr_Ccc.Tar 
               * Chr_Ccc.Apr + Chr_Ccc.Tpr * Chr_Ccc.Ppr );
 R_Cccc.Tpra = _d_Chr_Cccc.Tpar - _d_Chr_Cccc.Tpra + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpa + Chr_Ccc.Trr * Chr_Ccc.Rpa + Chr_Ccc.Tar 
               * Chr_Ccc.Apa + Chr_Ccc.Tpr * Chr_Ccc.Ppa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpr + Chr_Ccc.Tra * Chr_Ccc.Rpr + Chr_Ccc.Taa 
               * Chr_Ccc.Apr + Chr_Ccc.Tpa * Chr_Ccc.Ppr );
 R_Cccc.Tprp = _d_Chr_Cccc.Tppr - _d_Chr_Cccc.Tprp + ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpp + Chr_Ccc.Trr * Chr_Ccc.Rpp + Chr_Ccc.Tar 
               * Chr_Ccc.App + Chr_Ccc.Tpr * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpr + Chr_Ccc.Trp * Chr_Ccc.Rpr + Chr_Ccc.Tap 
               * Chr_Ccc.Apr + Chr_Ccc.Tpp * Chr_Ccc.Ppr );
 R_Cccc.Tpat = _d_Chr_Cccc.Tpta - _d_Chr_Cccc.Tpat + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpt + Chr_Ccc.Tra * Chr_Ccc.Rpt + Chr_Ccc.Taa 
               * Chr_Ccc.Apt + Chr_Ccc.Tpa * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpa + Chr_Ccc.Trt * Chr_Ccc.Rpa + Chr_Ccc.Tat 
               * Chr_Ccc.Apa + Chr_Ccc.Tpt * Chr_Ccc.Ppa );
 R_Cccc.Tpar = _d_Chr_Cccc.Tpra - _d_Chr_Cccc.Tpar + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpr + Chr_Ccc.Tra * Chr_Ccc.Rpr + Chr_Ccc.Taa 
               * Chr_Ccc.Apr + Chr_Ccc.Tpa * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpa + Chr_Ccc.Trr * Chr_Ccc.Rpa + Chr_Ccc.Tar 
               * Chr_Ccc.Apa + Chr_Ccc.Tpr * Chr_Ccc.Ppa );
 R_Cccc.Tpaa = _d_Chr_Cccc.Tpaa - _d_Chr_Cccc.Tpaa + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpa + Chr_Ccc.Tra * Chr_Ccc.Rpa + Chr_Ccc.Taa 
               * Chr_Ccc.Apa + Chr_Ccc.Tpa * Chr_Ccc.Ppa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpa + Chr_Ccc.Tra * Chr_Ccc.Rpa + Chr_Ccc.Taa 
               * Chr_Ccc.Apa + Chr_Ccc.Tpa * Chr_Ccc.Ppa );
 R_Cccc.Tpap = _d_Chr_Cccc.Tppa - _d_Chr_Cccc.Tpap + ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpp + Chr_Ccc.Tra * Chr_Ccc.Rpp + Chr_Ccc.Taa 
               * Chr_Ccc.App + Chr_Ccc.Tpa * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpa + Chr_Ccc.Trp * Chr_Ccc.Rpa + Chr_Ccc.Tap 
               * Chr_Ccc.Apa + Chr_Ccc.Tpp * Chr_Ccc.Ppa );
 R_Cccc.Tppt = _d_Chr_Cccc.Tptp - _d_Chr_Cccc.Tppt + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpt + Chr_Ccc.Trp * Chr_Ccc.Rpt + Chr_Ccc.Tap 
               * Chr_Ccc.Apt + Chr_Ccc.Tpp * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ttt 
               * Chr_Ccc.Tpp + Chr_Ccc.Trt * Chr_Ccc.Rpp + Chr_Ccc.Tat 
               * Chr_Ccc.App + Chr_Ccc.Tpt * Chr_Ccc.Ppp );
 R_Cccc.Tppr = _d_Chr_Cccc.Tprp - _d_Chr_Cccc.Tppr + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpr + Chr_Ccc.Trp * Chr_Ccc.Rpr + Chr_Ccc.Tap 
               * Chr_Ccc.Apr + Chr_Ccc.Tpp * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ttr 
               * Chr_Ccc.Tpp + Chr_Ccc.Trr * Chr_Ccc.Rpp + Chr_Ccc.Tar 
               * Chr_Ccc.App + Chr_Ccc.Tpr * Chr_Ccc.Ppp );
 R_Cccc.Tppa = _d_Chr_Cccc.Tpap - _d_Chr_Cccc.Tppa + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpa + Chr_Ccc.Trp * Chr_Ccc.Rpa + Chr_Ccc.Tap 
               * Chr_Ccc.Apa + Chr_Ccc.Tpp * Chr_Ccc.Ppa ) - ( Chr_Ccc.Tta 
               * Chr_Ccc.Tpp + Chr_Ccc.Tra * Chr_Ccc.Rpp + Chr_Ccc.Taa 
               * Chr_Ccc.App + Chr_Ccc.Tpa * Chr_Ccc.Ppp );
 R_Cccc.Tppp = _d_Chr_Cccc.Tppp - _d_Chr_Cccc.Tppp + ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpp + Chr_Ccc.Trp * Chr_Ccc.Rpp + Chr_Ccc.Tap 
               * Chr_Ccc.App + Chr_Ccc.Tpp * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ttp 
               * Chr_Ccc.Tpp + Chr_Ccc.Trp * Chr_Ccc.Rpp + Chr_Ccc.Tap 
               * Chr_Ccc.App + Chr_Ccc.Tpp * Chr_Ccc.Ppp );
 R_Cccc.Rttt = _d_Chr_Cccc.Rttt - _d_Chr_Cccc.Rttt + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Ttt + Chr_Ccc.Rrt * Chr_Ccc.Rtt + Chr_Ccc.Rat 
               * Chr_Ccc.Att + Chr_Ccc.Rpt * Chr_Ccc.Ptt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Ttt + Chr_Ccc.Rrt * Chr_Ccc.Rtt + Chr_Ccc.Rat 
               * Chr_Ccc.Att + Chr_Ccc.Rpt * Chr_Ccc.Ptt );
 R_Cccc.Rttr = _d_Chr_Cccc.Rtrt - _d_Chr_Cccc.Rttr + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Ttr + Chr_Ccc.Rrt * Chr_Ccc.Rtr + Chr_Ccc.Rat 
               * Chr_Ccc.Atr + Chr_Ccc.Rpt * Chr_Ccc.Ptr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Ttt + Chr_Ccc.Rrr * Chr_Ccc.Rtt + Chr_Ccc.Rar 
               * Chr_Ccc.Att + Chr_Ccc.Rpr * Chr_Ccc.Ptt );
 R_Cccc.Rtta = _d_Chr_Cccc.Rtat - _d_Chr_Cccc.Rtta + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tta + Chr_Ccc.Rrt * Chr_Ccc.Rta + Chr_Ccc.Rat 
               * Chr_Ccc.Ata + Chr_Ccc.Rpt * Chr_Ccc.Pta ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Ttt + Chr_Ccc.Rra * Chr_Ccc.Rtt + Chr_Ccc.Raa 
               * Chr_Ccc.Att + Chr_Ccc.Rpa * Chr_Ccc.Ptt );
 R_Cccc.Rttp = _d_Chr_Cccc.Rtpt - _d_Chr_Cccc.Rttp + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Ttp + Chr_Ccc.Rrt * Chr_Ccc.Rtp + Chr_Ccc.Rat 
               * Chr_Ccc.Atp + Chr_Ccc.Rpt * Chr_Ccc.Ptp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Ttt + Chr_Ccc.Rrp * Chr_Ccc.Rtt + Chr_Ccc.Rap 
               * Chr_Ccc.Att + Chr_Ccc.Rpp * Chr_Ccc.Ptt );
 R_Cccc.Rtrt = _d_Chr_Cccc.Rttr - _d_Chr_Cccc.Rtrt + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Ttt + Chr_Ccc.Rrr * Chr_Ccc.Rtt + Chr_Ccc.Rar 
               * Chr_Ccc.Att + Chr_Ccc.Rpr * Chr_Ccc.Ptt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Ttr + Chr_Ccc.Rrt * Chr_Ccc.Rtr + Chr_Ccc.Rat 
               * Chr_Ccc.Atr + Chr_Ccc.Rpt * Chr_Ccc.Ptr );
 R_Cccc.Rtrr = _d_Chr_Cccc.Rtrr - _d_Chr_Cccc.Rtrr + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Ttr + Chr_Ccc.Rrr * Chr_Ccc.Rtr + Chr_Ccc.Rar 
               * Chr_Ccc.Atr + Chr_Ccc.Rpr * Chr_Ccc.Ptr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Ttr + Chr_Ccc.Rrr * Chr_Ccc.Rtr + Chr_Ccc.Rar 
               * Chr_Ccc.Atr + Chr_Ccc.Rpr * Chr_Ccc.Ptr );
 R_Cccc.Rtra = _d_Chr_Cccc.Rtar - _d_Chr_Cccc.Rtra + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tta + Chr_Ccc.Rrr * Chr_Ccc.Rta + Chr_Ccc.Rar 
               * Chr_Ccc.Ata + Chr_Ccc.Rpr * Chr_Ccc.Pta ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Ttr + Chr_Ccc.Rra * Chr_Ccc.Rtr + Chr_Ccc.Raa 
               * Chr_Ccc.Atr + Chr_Ccc.Rpa * Chr_Ccc.Ptr );
 R_Cccc.Rtrp = _d_Chr_Cccc.Rtpr - _d_Chr_Cccc.Rtrp + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Ttp + Chr_Ccc.Rrr * Chr_Ccc.Rtp + Chr_Ccc.Rar 
               * Chr_Ccc.Atp + Chr_Ccc.Rpr * Chr_Ccc.Ptp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Ttr + Chr_Ccc.Rrp * Chr_Ccc.Rtr + Chr_Ccc.Rap 
               * Chr_Ccc.Atr + Chr_Ccc.Rpp * Chr_Ccc.Ptr );
 R_Cccc.Rtat = _d_Chr_Cccc.Rtta - _d_Chr_Cccc.Rtat + ( Chr_Ccc.Rta 
               * Chr_Ccc.Ttt + Chr_Ccc.Rra * Chr_Ccc.Rtt + Chr_Ccc.Raa 
               * Chr_Ccc.Att + Chr_Ccc.Rpa * Chr_Ccc.Ptt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tta + Chr_Ccc.Rrt * Chr_Ccc.Rta + Chr_Ccc.Rat 
               * Chr_Ccc.Ata + Chr_Ccc.Rpt * Chr_Ccc.Pta );
 R_Cccc.Rtar = _d_Chr_Cccc.Rtra - _d_Chr_Cccc.Rtar + ( Chr_Ccc.Rta 
               * Chr_Ccc.Ttr + Chr_Ccc.Rra * Chr_Ccc.Rtr + Chr_Ccc.Raa 
               * Chr_Ccc.Atr + Chr_Ccc.Rpa * Chr_Ccc.Ptr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tta + Chr_Ccc.Rrr * Chr_Ccc.Rta + Chr_Ccc.Rar 
               * Chr_Ccc.Ata + Chr_Ccc.Rpr * Chr_Ccc.Pta );
 R_Cccc.Rtaa = _d_Chr_Cccc.Rtaa - _d_Chr_Cccc.Rtaa + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tta + Chr_Ccc.Rra * Chr_Ccc.Rta + Chr_Ccc.Raa 
               * Chr_Ccc.Ata + Chr_Ccc.Rpa * Chr_Ccc.Pta ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tta + Chr_Ccc.Rra * Chr_Ccc.Rta + Chr_Ccc.Raa 
               * Chr_Ccc.Ata + Chr_Ccc.Rpa * Chr_Ccc.Pta );
 R_Cccc.Rtap = _d_Chr_Cccc.Rtpa - _d_Chr_Cccc.Rtap + ( Chr_Ccc.Rta 
               * Chr_Ccc.Ttp + Chr_Ccc.Rra * Chr_Ccc.Rtp + Chr_Ccc.Raa 
               * Chr_Ccc.Atp + Chr_Ccc.Rpa * Chr_Ccc.Ptp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tta + Chr_Ccc.Rrp * Chr_Ccc.Rta + Chr_Ccc.Rap 
               * Chr_Ccc.Ata + Chr_Ccc.Rpp * Chr_Ccc.Pta );
 R_Cccc.Rtpt = _d_Chr_Cccc.Rttp - _d_Chr_Cccc.Rtpt + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Ttt + Chr_Ccc.Rrp * Chr_Ccc.Rtt + Chr_Ccc.Rap 
               * Chr_Ccc.Att + Chr_Ccc.Rpp * Chr_Ccc.Ptt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Ttp + Chr_Ccc.Rrt * Chr_Ccc.Rtp + Chr_Ccc.Rat 
               * Chr_Ccc.Atp + Chr_Ccc.Rpt * Chr_Ccc.Ptp );
 R_Cccc.Rtpr = _d_Chr_Cccc.Rtrp - _d_Chr_Cccc.Rtpr + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Ttr + Chr_Ccc.Rrp * Chr_Ccc.Rtr + Chr_Ccc.Rap 
               * Chr_Ccc.Atr + Chr_Ccc.Rpp * Chr_Ccc.Ptr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Ttp + Chr_Ccc.Rrr * Chr_Ccc.Rtp + Chr_Ccc.Rar 
               * Chr_Ccc.Atp + Chr_Ccc.Rpr * Chr_Ccc.Ptp );
 R_Cccc.Rtpa = _d_Chr_Cccc.Rtap - _d_Chr_Cccc.Rtpa + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tta + Chr_Ccc.Rrp * Chr_Ccc.Rta + Chr_Ccc.Rap 
               * Chr_Ccc.Ata + Chr_Ccc.Rpp * Chr_Ccc.Pta ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Ttp + Chr_Ccc.Rra * Chr_Ccc.Rtp + Chr_Ccc.Raa 
               * Chr_Ccc.Atp + Chr_Ccc.Rpa * Chr_Ccc.Ptp );
 R_Cccc.Rtpp = _d_Chr_Cccc.Rtpp - _d_Chr_Cccc.Rtpp + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Ttp + Chr_Ccc.Rrp * Chr_Ccc.Rtp + Chr_Ccc.Rap 
               * Chr_Ccc.Atp + Chr_Ccc.Rpp * Chr_Ccc.Ptp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Ttp + Chr_Ccc.Rrp * Chr_Ccc.Rtp + Chr_Ccc.Rap 
               * Chr_Ccc.Atp + Chr_Ccc.Rpp * Chr_Ccc.Ptp );
 R_Cccc.Rrtt = _d_Chr_Cccc.Rrtt - _d_Chr_Cccc.Rrtt + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Trt + Chr_Ccc.Rrt * Chr_Ccc.Rrt + Chr_Ccc.Rat 
               * Chr_Ccc.Art + Chr_Ccc.Rpt * Chr_Ccc.Prt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Trt + Chr_Ccc.Rrt * Chr_Ccc.Rrt + Chr_Ccc.Rat 
               * Chr_Ccc.Art + Chr_Ccc.Rpt * Chr_Ccc.Prt );
 R_Cccc.Rrtr = _d_Chr_Cccc.Rrrt - _d_Chr_Cccc.Rrtr + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Trr + Chr_Ccc.Rrt * Chr_Ccc.Rrr + Chr_Ccc.Rat 
               * Chr_Ccc.Arr + Chr_Ccc.Rpt * Chr_Ccc.Prr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Trt + Chr_Ccc.Rrr * Chr_Ccc.Rrt + Chr_Ccc.Rar 
               * Chr_Ccc.Art + Chr_Ccc.Rpr * Chr_Ccc.Prt );
 R_Cccc.Rrta = _d_Chr_Cccc.Rrat - _d_Chr_Cccc.Rrta + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tra + Chr_Ccc.Rrt * Chr_Ccc.Rra + Chr_Ccc.Rat 
               * Chr_Ccc.Ara + Chr_Ccc.Rpt * Chr_Ccc.Pra ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Trt + Chr_Ccc.Rra * Chr_Ccc.Rrt + Chr_Ccc.Raa 
               * Chr_Ccc.Art + Chr_Ccc.Rpa * Chr_Ccc.Prt );
 R_Cccc.Rrtp = _d_Chr_Cccc.Rrpt - _d_Chr_Cccc.Rrtp + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Trp + Chr_Ccc.Rrt * Chr_Ccc.Rrp + Chr_Ccc.Rat 
               * Chr_Ccc.Arp + Chr_Ccc.Rpt * Chr_Ccc.Prp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Trt + Chr_Ccc.Rrp * Chr_Ccc.Rrt + Chr_Ccc.Rap 
               * Chr_Ccc.Art + Chr_Ccc.Rpp * Chr_Ccc.Prt );
 R_Cccc.Rrrt = _d_Chr_Cccc.Rrtr - _d_Chr_Cccc.Rrrt + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Trt + Chr_Ccc.Rrr * Chr_Ccc.Rrt + Chr_Ccc.Rar 
               * Chr_Ccc.Art + Chr_Ccc.Rpr * Chr_Ccc.Prt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Trr + Chr_Ccc.Rrt * Chr_Ccc.Rrr + Chr_Ccc.Rat 
               * Chr_Ccc.Arr + Chr_Ccc.Rpt * Chr_Ccc.Prr );
 R_Cccc.Rrrr = _d_Chr_Cccc.Rrrr - _d_Chr_Cccc.Rrrr + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Trr + Chr_Ccc.Rrr * Chr_Ccc.Rrr + Chr_Ccc.Rar 
               * Chr_Ccc.Arr + Chr_Ccc.Rpr * Chr_Ccc.Prr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Trr + Chr_Ccc.Rrr * Chr_Ccc.Rrr + Chr_Ccc.Rar 
               * Chr_Ccc.Arr + Chr_Ccc.Rpr * Chr_Ccc.Prr );
 R_Cccc.Rrra = _d_Chr_Cccc.Rrar - _d_Chr_Cccc.Rrra + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tra + Chr_Ccc.Rrr * Chr_Ccc.Rra + Chr_Ccc.Rar 
               * Chr_Ccc.Ara + Chr_Ccc.Rpr * Chr_Ccc.Pra ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Trr + Chr_Ccc.Rra * Chr_Ccc.Rrr + Chr_Ccc.Raa 
               * Chr_Ccc.Arr + Chr_Ccc.Rpa * Chr_Ccc.Prr );
 R_Cccc.Rrrp = _d_Chr_Cccc.Rrpr - _d_Chr_Cccc.Rrrp + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Trp + Chr_Ccc.Rrr * Chr_Ccc.Rrp + Chr_Ccc.Rar 
               * Chr_Ccc.Arp + Chr_Ccc.Rpr * Chr_Ccc.Prp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Trr + Chr_Ccc.Rrp * Chr_Ccc.Rrr + Chr_Ccc.Rap 
               * Chr_Ccc.Arr + Chr_Ccc.Rpp * Chr_Ccc.Prr );
 R_Cccc.Rrat = _d_Chr_Cccc.Rrta - _d_Chr_Cccc.Rrat + ( Chr_Ccc.Rta 
               * Chr_Ccc.Trt + Chr_Ccc.Rra * Chr_Ccc.Rrt + Chr_Ccc.Raa 
               * Chr_Ccc.Art + Chr_Ccc.Rpa * Chr_Ccc.Prt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tra + Chr_Ccc.Rrt * Chr_Ccc.Rra + Chr_Ccc.Rat 
               * Chr_Ccc.Ara + Chr_Ccc.Rpt * Chr_Ccc.Pra );
 R_Cccc.Rrar = _d_Chr_Cccc.Rrra - _d_Chr_Cccc.Rrar + ( Chr_Ccc.Rta 
               * Chr_Ccc.Trr + Chr_Ccc.Rra * Chr_Ccc.Rrr + Chr_Ccc.Raa 
               * Chr_Ccc.Arr + Chr_Ccc.Rpa * Chr_Ccc.Prr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tra + Chr_Ccc.Rrr * Chr_Ccc.Rra + Chr_Ccc.Rar 
               * Chr_Ccc.Ara + Chr_Ccc.Rpr * Chr_Ccc.Pra );
 R_Cccc.Rraa = _d_Chr_Cccc.Rraa - _d_Chr_Cccc.Rraa + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tra + Chr_Ccc.Rra * Chr_Ccc.Rra + Chr_Ccc.Raa 
               * Chr_Ccc.Ara + Chr_Ccc.Rpa * Chr_Ccc.Pra ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tra + Chr_Ccc.Rra * Chr_Ccc.Rra + Chr_Ccc.Raa 
               * Chr_Ccc.Ara + Chr_Ccc.Rpa * Chr_Ccc.Pra );
 R_Cccc.Rrap = _d_Chr_Cccc.Rrpa - _d_Chr_Cccc.Rrap + ( Chr_Ccc.Rta 
               * Chr_Ccc.Trp + Chr_Ccc.Rra * Chr_Ccc.Rrp + Chr_Ccc.Raa 
               * Chr_Ccc.Arp + Chr_Ccc.Rpa * Chr_Ccc.Prp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tra + Chr_Ccc.Rrp * Chr_Ccc.Rra + Chr_Ccc.Rap 
               * Chr_Ccc.Ara + Chr_Ccc.Rpp * Chr_Ccc.Pra );
 R_Cccc.Rrpt = _d_Chr_Cccc.Rrtp - _d_Chr_Cccc.Rrpt + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Trt + Chr_Ccc.Rrp * Chr_Ccc.Rrt + Chr_Ccc.Rap 
               * Chr_Ccc.Art + Chr_Ccc.Rpp * Chr_Ccc.Prt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Trp + Chr_Ccc.Rrt * Chr_Ccc.Rrp + Chr_Ccc.Rat 
               * Chr_Ccc.Arp + Chr_Ccc.Rpt * Chr_Ccc.Prp );
 R_Cccc.Rrpr = _d_Chr_Cccc.Rrrp - _d_Chr_Cccc.Rrpr + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Trr + Chr_Ccc.Rrp * Chr_Ccc.Rrr + Chr_Ccc.Rap 
               * Chr_Ccc.Arr + Chr_Ccc.Rpp * Chr_Ccc.Prr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Trp + Chr_Ccc.Rrr * Chr_Ccc.Rrp + Chr_Ccc.Rar 
               * Chr_Ccc.Arp + Chr_Ccc.Rpr * Chr_Ccc.Prp );
 R_Cccc.Rrpa = _d_Chr_Cccc.Rrap - _d_Chr_Cccc.Rrpa + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tra + Chr_Ccc.Rrp * Chr_Ccc.Rra + Chr_Ccc.Rap 
               * Chr_Ccc.Ara + Chr_Ccc.Rpp * Chr_Ccc.Pra ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Trp + Chr_Ccc.Rra * Chr_Ccc.Rrp + Chr_Ccc.Raa 
               * Chr_Ccc.Arp + Chr_Ccc.Rpa * Chr_Ccc.Prp );
 R_Cccc.Rrpp = _d_Chr_Cccc.Rrpp - _d_Chr_Cccc.Rrpp + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Trp + Chr_Ccc.Rrp * Chr_Ccc.Rrp + Chr_Ccc.Rap 
               * Chr_Ccc.Arp + Chr_Ccc.Rpp * Chr_Ccc.Prp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Trp + Chr_Ccc.Rrp * Chr_Ccc.Rrp + Chr_Ccc.Rap 
               * Chr_Ccc.Arp + Chr_Ccc.Rpp * Chr_Ccc.Prp );
 R_Cccc.Ratt = _d_Chr_Cccc.Ratt - _d_Chr_Cccc.Ratt + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tat + Chr_Ccc.Rrt * Chr_Ccc.Rat + Chr_Ccc.Rat 
               * Chr_Ccc.Aat + Chr_Ccc.Rpt * Chr_Ccc.Pat ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tat + Chr_Ccc.Rrt * Chr_Ccc.Rat + Chr_Ccc.Rat 
               * Chr_Ccc.Aat + Chr_Ccc.Rpt * Chr_Ccc.Pat );
 R_Cccc.Ratr = _d_Chr_Cccc.Rart - _d_Chr_Cccc.Ratr + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tar + Chr_Ccc.Rrt * Chr_Ccc.Rar + Chr_Ccc.Rat 
               * Chr_Ccc.Aar + Chr_Ccc.Rpt * Chr_Ccc.Par ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tat + Chr_Ccc.Rrr * Chr_Ccc.Rat + Chr_Ccc.Rar 
               * Chr_Ccc.Aat + Chr_Ccc.Rpr * Chr_Ccc.Pat );
 R_Cccc.Rata = _d_Chr_Cccc.Raat - _d_Chr_Cccc.Rata + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Taa + Chr_Ccc.Rrt * Chr_Ccc.Raa + Chr_Ccc.Rat 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpt * Chr_Ccc.Paa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tat + Chr_Ccc.Rra * Chr_Ccc.Rat + Chr_Ccc.Raa 
               * Chr_Ccc.Aat + Chr_Ccc.Rpa * Chr_Ccc.Pat );
 R_Cccc.Ratp = _d_Chr_Cccc.Rapt - _d_Chr_Cccc.Ratp + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tap + Chr_Ccc.Rrt * Chr_Ccc.Rap + Chr_Ccc.Rat 
               * Chr_Ccc.Aap + Chr_Ccc.Rpt * Chr_Ccc.Pap ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tat + Chr_Ccc.Rrp * Chr_Ccc.Rat + Chr_Ccc.Rap 
               * Chr_Ccc.Aat + Chr_Ccc.Rpp * Chr_Ccc.Pat );
 R_Cccc.Rart = _d_Chr_Cccc.Ratr - _d_Chr_Cccc.Rart + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tat + Chr_Ccc.Rrr * Chr_Ccc.Rat + Chr_Ccc.Rar 
               * Chr_Ccc.Aat + Chr_Ccc.Rpr * Chr_Ccc.Pat ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tar + Chr_Ccc.Rrt * Chr_Ccc.Rar + Chr_Ccc.Rat 
               * Chr_Ccc.Aar + Chr_Ccc.Rpt * Chr_Ccc.Par );
 R_Cccc.Rarr = _d_Chr_Cccc.Rarr - _d_Chr_Cccc.Rarr + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tar + Chr_Ccc.Rrr * Chr_Ccc.Rar + Chr_Ccc.Rar 
               * Chr_Ccc.Aar + Chr_Ccc.Rpr * Chr_Ccc.Par ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tar + Chr_Ccc.Rrr * Chr_Ccc.Rar + Chr_Ccc.Rar 
               * Chr_Ccc.Aar + Chr_Ccc.Rpr * Chr_Ccc.Par );
 R_Cccc.Rara = _d_Chr_Cccc.Raar - _d_Chr_Cccc.Rara + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Taa + Chr_Ccc.Rrr * Chr_Ccc.Raa + Chr_Ccc.Rar 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpr * Chr_Ccc.Paa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tar + Chr_Ccc.Rra * Chr_Ccc.Rar + Chr_Ccc.Raa 
               * Chr_Ccc.Aar + Chr_Ccc.Rpa * Chr_Ccc.Par );
 R_Cccc.Rarp = _d_Chr_Cccc.Rapr - _d_Chr_Cccc.Rarp + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tap + Chr_Ccc.Rrr * Chr_Ccc.Rap + Chr_Ccc.Rar 
               * Chr_Ccc.Aap + Chr_Ccc.Rpr * Chr_Ccc.Pap ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tar + Chr_Ccc.Rrp * Chr_Ccc.Rar + Chr_Ccc.Rap 
               * Chr_Ccc.Aar + Chr_Ccc.Rpp * Chr_Ccc.Par );
 R_Cccc.Raat = _d_Chr_Cccc.Rata - _d_Chr_Cccc.Raat + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tat + Chr_Ccc.Rra * Chr_Ccc.Rat + Chr_Ccc.Raa 
               * Chr_Ccc.Aat + Chr_Ccc.Rpa * Chr_Ccc.Pat ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Taa + Chr_Ccc.Rrt * Chr_Ccc.Raa + Chr_Ccc.Rat 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpt * Chr_Ccc.Paa );
 R_Cccc.Raar = _d_Chr_Cccc.Rara - _d_Chr_Cccc.Raar + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tar + Chr_Ccc.Rra * Chr_Ccc.Rar + Chr_Ccc.Raa 
               * Chr_Ccc.Aar + Chr_Ccc.Rpa * Chr_Ccc.Par ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Taa + Chr_Ccc.Rrr * Chr_Ccc.Raa + Chr_Ccc.Rar 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpr * Chr_Ccc.Paa );
 R_Cccc.Raaa = _d_Chr_Cccc.Raaa - _d_Chr_Cccc.Raaa + ( Chr_Ccc.Rta 
               * Chr_Ccc.Taa + Chr_Ccc.Rra * Chr_Ccc.Raa + Chr_Ccc.Raa 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpa * Chr_Ccc.Paa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Taa + Chr_Ccc.Rra * Chr_Ccc.Raa + Chr_Ccc.Raa 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpa * Chr_Ccc.Paa );
 R_Cccc.Raap = _d_Chr_Cccc.Rapa - _d_Chr_Cccc.Raap + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tap + Chr_Ccc.Rra * Chr_Ccc.Rap + Chr_Ccc.Raa 
               * Chr_Ccc.Aap + Chr_Ccc.Rpa * Chr_Ccc.Pap ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Taa + Chr_Ccc.Rrp * Chr_Ccc.Raa + Chr_Ccc.Rap 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpp * Chr_Ccc.Paa );
 R_Cccc.Rapt = _d_Chr_Cccc.Ratp - _d_Chr_Cccc.Rapt + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tat + Chr_Ccc.Rrp * Chr_Ccc.Rat + Chr_Ccc.Rap 
               * Chr_Ccc.Aat + Chr_Ccc.Rpp * Chr_Ccc.Pat ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tap + Chr_Ccc.Rrt * Chr_Ccc.Rap + Chr_Ccc.Rat 
               * Chr_Ccc.Aap + Chr_Ccc.Rpt * Chr_Ccc.Pap );
 R_Cccc.Rapr = _d_Chr_Cccc.Rarp - _d_Chr_Cccc.Rapr + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tar + Chr_Ccc.Rrp * Chr_Ccc.Rar + Chr_Ccc.Rap 
               * Chr_Ccc.Aar + Chr_Ccc.Rpp * Chr_Ccc.Par ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tap + Chr_Ccc.Rrr * Chr_Ccc.Rap + Chr_Ccc.Rar 
               * Chr_Ccc.Aap + Chr_Ccc.Rpr * Chr_Ccc.Pap );
 R_Cccc.Rapa = _d_Chr_Cccc.Raap - _d_Chr_Cccc.Rapa + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Taa + Chr_Ccc.Rrp * Chr_Ccc.Raa + Chr_Ccc.Rap 
               * Chr_Ccc.Aaa + Chr_Ccc.Rpp * Chr_Ccc.Paa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tap + Chr_Ccc.Rra * Chr_Ccc.Rap + Chr_Ccc.Raa 
               * Chr_Ccc.Aap + Chr_Ccc.Rpa * Chr_Ccc.Pap );
 R_Cccc.Rapp = _d_Chr_Cccc.Rapp - _d_Chr_Cccc.Rapp + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tap + Chr_Ccc.Rrp * Chr_Ccc.Rap + Chr_Ccc.Rap 
               * Chr_Ccc.Aap + Chr_Ccc.Rpp * Chr_Ccc.Pap ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tap + Chr_Ccc.Rrp * Chr_Ccc.Rap + Chr_Ccc.Rap 
               * Chr_Ccc.Aap + Chr_Ccc.Rpp * Chr_Ccc.Pap );
 R_Cccc.Rptt = _d_Chr_Cccc.Rptt - _d_Chr_Cccc.Rptt + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpt + Chr_Ccc.Rrt * Chr_Ccc.Rpt + Chr_Ccc.Rat 
               * Chr_Ccc.Apt + Chr_Ccc.Rpt * Chr_Ccc.Ppt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpt + Chr_Ccc.Rrt * Chr_Ccc.Rpt + Chr_Ccc.Rat 
               * Chr_Ccc.Apt + Chr_Ccc.Rpt * Chr_Ccc.Ppt );
 R_Cccc.Rptr = _d_Chr_Cccc.Rprt - _d_Chr_Cccc.Rptr + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpr + Chr_Ccc.Rrt * Chr_Ccc.Rpr + Chr_Ccc.Rat 
               * Chr_Ccc.Apr + Chr_Ccc.Rpt * Chr_Ccc.Ppr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpt + Chr_Ccc.Rrr * Chr_Ccc.Rpt + Chr_Ccc.Rar 
               * Chr_Ccc.Apt + Chr_Ccc.Rpr * Chr_Ccc.Ppt );
 R_Cccc.Rpta = _d_Chr_Cccc.Rpat - _d_Chr_Cccc.Rpta + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpa + Chr_Ccc.Rrt * Chr_Ccc.Rpa + Chr_Ccc.Rat 
               * Chr_Ccc.Apa + Chr_Ccc.Rpt * Chr_Ccc.Ppa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpt + Chr_Ccc.Rra * Chr_Ccc.Rpt + Chr_Ccc.Raa 
               * Chr_Ccc.Apt + Chr_Ccc.Rpa * Chr_Ccc.Ppt );
 R_Cccc.Rptp = _d_Chr_Cccc.Rppt - _d_Chr_Cccc.Rptp + ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpp + Chr_Ccc.Rrt * Chr_Ccc.Rpp + Chr_Ccc.Rat 
               * Chr_Ccc.App + Chr_Ccc.Rpt * Chr_Ccc.Ppp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpt + Chr_Ccc.Rrp * Chr_Ccc.Rpt + Chr_Ccc.Rap 
               * Chr_Ccc.Apt + Chr_Ccc.Rpp * Chr_Ccc.Ppt );
 R_Cccc.Rprt = _d_Chr_Cccc.Rptr - _d_Chr_Cccc.Rprt + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpt + Chr_Ccc.Rrr * Chr_Ccc.Rpt + Chr_Ccc.Rar 
               * Chr_Ccc.Apt + Chr_Ccc.Rpr * Chr_Ccc.Ppt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpr + Chr_Ccc.Rrt * Chr_Ccc.Rpr + Chr_Ccc.Rat 
               * Chr_Ccc.Apr + Chr_Ccc.Rpt * Chr_Ccc.Ppr );
 R_Cccc.Rprr = _d_Chr_Cccc.Rprr - _d_Chr_Cccc.Rprr + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpr + Chr_Ccc.Rrr * Chr_Ccc.Rpr + Chr_Ccc.Rar 
               * Chr_Ccc.Apr + Chr_Ccc.Rpr * Chr_Ccc.Ppr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpr + Chr_Ccc.Rrr * Chr_Ccc.Rpr + Chr_Ccc.Rar 
               * Chr_Ccc.Apr + Chr_Ccc.Rpr * Chr_Ccc.Ppr );
 R_Cccc.Rpra = _d_Chr_Cccc.Rpar - _d_Chr_Cccc.Rpra + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpa + Chr_Ccc.Rrr * Chr_Ccc.Rpa + Chr_Ccc.Rar 
               * Chr_Ccc.Apa + Chr_Ccc.Rpr * Chr_Ccc.Ppa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpr + Chr_Ccc.Rra * Chr_Ccc.Rpr + Chr_Ccc.Raa 
               * Chr_Ccc.Apr + Chr_Ccc.Rpa * Chr_Ccc.Ppr );
 R_Cccc.Rprp = _d_Chr_Cccc.Rppr - _d_Chr_Cccc.Rprp + ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpp + Chr_Ccc.Rrr * Chr_Ccc.Rpp + Chr_Ccc.Rar 
               * Chr_Ccc.App + Chr_Ccc.Rpr * Chr_Ccc.Ppp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpr + Chr_Ccc.Rrp * Chr_Ccc.Rpr + Chr_Ccc.Rap 
               * Chr_Ccc.Apr + Chr_Ccc.Rpp * Chr_Ccc.Ppr );
 R_Cccc.Rpat = _d_Chr_Cccc.Rpta - _d_Chr_Cccc.Rpat + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpt + Chr_Ccc.Rra * Chr_Ccc.Rpt + Chr_Ccc.Raa 
               * Chr_Ccc.Apt + Chr_Ccc.Rpa * Chr_Ccc.Ppt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpa + Chr_Ccc.Rrt * Chr_Ccc.Rpa + Chr_Ccc.Rat 
               * Chr_Ccc.Apa + Chr_Ccc.Rpt * Chr_Ccc.Ppa );
 R_Cccc.Rpar = _d_Chr_Cccc.Rpra - _d_Chr_Cccc.Rpar + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpr + Chr_Ccc.Rra * Chr_Ccc.Rpr + Chr_Ccc.Raa 
               * Chr_Ccc.Apr + Chr_Ccc.Rpa * Chr_Ccc.Ppr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpa + Chr_Ccc.Rrr * Chr_Ccc.Rpa + Chr_Ccc.Rar 
               * Chr_Ccc.Apa + Chr_Ccc.Rpr * Chr_Ccc.Ppa );
 R_Cccc.Rpaa = _d_Chr_Cccc.Rpaa - _d_Chr_Cccc.Rpaa + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpa + Chr_Ccc.Rra * Chr_Ccc.Rpa + Chr_Ccc.Raa 
               * Chr_Ccc.Apa + Chr_Ccc.Rpa * Chr_Ccc.Ppa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpa + Chr_Ccc.Rra * Chr_Ccc.Rpa + Chr_Ccc.Raa 
               * Chr_Ccc.Apa + Chr_Ccc.Rpa * Chr_Ccc.Ppa );
 R_Cccc.Rpap = _d_Chr_Cccc.Rppa - _d_Chr_Cccc.Rpap + ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpp + Chr_Ccc.Rra * Chr_Ccc.Rpp + Chr_Ccc.Raa 
               * Chr_Ccc.App + Chr_Ccc.Rpa * Chr_Ccc.Ppp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpa + Chr_Ccc.Rrp * Chr_Ccc.Rpa + Chr_Ccc.Rap 
               * Chr_Ccc.Apa + Chr_Ccc.Rpp * Chr_Ccc.Ppa );
 R_Cccc.Rppt = _d_Chr_Cccc.Rptp - _d_Chr_Cccc.Rppt + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpt + Chr_Ccc.Rrp * Chr_Ccc.Rpt + Chr_Ccc.Rap 
               * Chr_Ccc.Apt + Chr_Ccc.Rpp * Chr_Ccc.Ppt ) - ( Chr_Ccc.Rtt 
               * Chr_Ccc.Tpp + Chr_Ccc.Rrt * Chr_Ccc.Rpp + Chr_Ccc.Rat 
               * Chr_Ccc.App + Chr_Ccc.Rpt * Chr_Ccc.Ppp );
 R_Cccc.Rppr = _d_Chr_Cccc.Rprp - _d_Chr_Cccc.Rppr + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpr + Chr_Ccc.Rrp * Chr_Ccc.Rpr + Chr_Ccc.Rap 
               * Chr_Ccc.Apr + Chr_Ccc.Rpp * Chr_Ccc.Ppr ) - ( Chr_Ccc.Rtr 
               * Chr_Ccc.Tpp + Chr_Ccc.Rrr * Chr_Ccc.Rpp + Chr_Ccc.Rar 
               * Chr_Ccc.App + Chr_Ccc.Rpr * Chr_Ccc.Ppp );
 R_Cccc.Rppa = _d_Chr_Cccc.Rpap - _d_Chr_Cccc.Rppa + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpa + Chr_Ccc.Rrp * Chr_Ccc.Rpa + Chr_Ccc.Rap 
               * Chr_Ccc.Apa + Chr_Ccc.Rpp * Chr_Ccc.Ppa ) - ( Chr_Ccc.Rta 
               * Chr_Ccc.Tpp + Chr_Ccc.Rra * Chr_Ccc.Rpp + Chr_Ccc.Raa 
               * Chr_Ccc.App + Chr_Ccc.Rpa * Chr_Ccc.Ppp );
 R_Cccc.Rppp = _d_Chr_Cccc.Rppp - _d_Chr_Cccc.Rppp + ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpp + Chr_Ccc.Rrp * Chr_Ccc.Rpp + Chr_Ccc.Rap 
               * Chr_Ccc.App + Chr_Ccc.Rpp * Chr_Ccc.Ppp ) - ( Chr_Ccc.Rtp 
               * Chr_Ccc.Tpp + Chr_Ccc.Rrp * Chr_Ccc.Rpp + Chr_Ccc.Rap 
               * Chr_Ccc.App + Chr_Ccc.Rpp * Chr_Ccc.Ppp );
 R_Cccc.Attt = _d_Chr_Cccc.Attt - _d_Chr_Cccc.Attt + ( Chr_Ccc.Att 
               * Chr_Ccc.Ttt + Chr_Ccc.Art * Chr_Ccc.Rtt + Chr_Ccc.Aat 
               * Chr_Ccc.Att + Chr_Ccc.Apt * Chr_Ccc.Ptt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Ttt + Chr_Ccc.Art * Chr_Ccc.Rtt + Chr_Ccc.Aat 
               * Chr_Ccc.Att + Chr_Ccc.Apt * Chr_Ccc.Ptt );
 R_Cccc.Attr = _d_Chr_Cccc.Atrt - _d_Chr_Cccc.Attr + ( Chr_Ccc.Att 
               * Chr_Ccc.Ttr + Chr_Ccc.Art * Chr_Ccc.Rtr + Chr_Ccc.Aat 
               * Chr_Ccc.Atr + Chr_Ccc.Apt * Chr_Ccc.Ptr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Ttt + Chr_Ccc.Arr * Chr_Ccc.Rtt + Chr_Ccc.Aar 
               * Chr_Ccc.Att + Chr_Ccc.Apr * Chr_Ccc.Ptt );
 R_Cccc.Atta = _d_Chr_Cccc.Atat - _d_Chr_Cccc.Atta + ( Chr_Ccc.Att 
               * Chr_Ccc.Tta + Chr_Ccc.Art * Chr_Ccc.Rta + Chr_Ccc.Aat 
               * Chr_Ccc.Ata + Chr_Ccc.Apt * Chr_Ccc.Pta ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Ttt + Chr_Ccc.Ara * Chr_Ccc.Rtt + Chr_Ccc.Aaa 
               * Chr_Ccc.Att + Chr_Ccc.Apa * Chr_Ccc.Ptt );
 R_Cccc.Attp = _d_Chr_Cccc.Atpt - _d_Chr_Cccc.Attp + ( Chr_Ccc.Att 
               * Chr_Ccc.Ttp + Chr_Ccc.Art * Chr_Ccc.Rtp + Chr_Ccc.Aat 
               * Chr_Ccc.Atp + Chr_Ccc.Apt * Chr_Ccc.Ptp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Ttt + Chr_Ccc.Arp * Chr_Ccc.Rtt + Chr_Ccc.Aap 
               * Chr_Ccc.Att + Chr_Ccc.App * Chr_Ccc.Ptt );
 R_Cccc.Atrt = _d_Chr_Cccc.Attr - _d_Chr_Cccc.Atrt + ( Chr_Ccc.Atr 
               * Chr_Ccc.Ttt + Chr_Ccc.Arr * Chr_Ccc.Rtt + Chr_Ccc.Aar 
               * Chr_Ccc.Att + Chr_Ccc.Apr * Chr_Ccc.Ptt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Ttr + Chr_Ccc.Art * Chr_Ccc.Rtr + Chr_Ccc.Aat 
               * Chr_Ccc.Atr + Chr_Ccc.Apt * Chr_Ccc.Ptr );
 R_Cccc.Atrr = _d_Chr_Cccc.Atrr - _d_Chr_Cccc.Atrr + ( Chr_Ccc.Atr 
               * Chr_Ccc.Ttr + Chr_Ccc.Arr * Chr_Ccc.Rtr + Chr_Ccc.Aar 
               * Chr_Ccc.Atr + Chr_Ccc.Apr * Chr_Ccc.Ptr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Ttr + Chr_Ccc.Arr * Chr_Ccc.Rtr + Chr_Ccc.Aar 
               * Chr_Ccc.Atr + Chr_Ccc.Apr * Chr_Ccc.Ptr );
 R_Cccc.Atra = _d_Chr_Cccc.Atar - _d_Chr_Cccc.Atra + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tta + Chr_Ccc.Arr * Chr_Ccc.Rta + Chr_Ccc.Aar 
               * Chr_Ccc.Ata + Chr_Ccc.Apr * Chr_Ccc.Pta ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Ttr + Chr_Ccc.Ara * Chr_Ccc.Rtr + Chr_Ccc.Aaa 
               * Chr_Ccc.Atr + Chr_Ccc.Apa * Chr_Ccc.Ptr );
 R_Cccc.Atrp = _d_Chr_Cccc.Atpr - _d_Chr_Cccc.Atrp + ( Chr_Ccc.Atr 
               * Chr_Ccc.Ttp + Chr_Ccc.Arr * Chr_Ccc.Rtp + Chr_Ccc.Aar 
               * Chr_Ccc.Atp + Chr_Ccc.Apr * Chr_Ccc.Ptp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Ttr + Chr_Ccc.Arp * Chr_Ccc.Rtr + Chr_Ccc.Aap 
               * Chr_Ccc.Atr + Chr_Ccc.App * Chr_Ccc.Ptr );
 R_Cccc.Atat = _d_Chr_Cccc.Atta - _d_Chr_Cccc.Atat + ( Chr_Ccc.Ata 
               * Chr_Ccc.Ttt + Chr_Ccc.Ara * Chr_Ccc.Rtt + Chr_Ccc.Aaa 
               * Chr_Ccc.Att + Chr_Ccc.Apa * Chr_Ccc.Ptt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tta + Chr_Ccc.Art * Chr_Ccc.Rta + Chr_Ccc.Aat 
               * Chr_Ccc.Ata + Chr_Ccc.Apt * Chr_Ccc.Pta );
 R_Cccc.Atar = _d_Chr_Cccc.Atra - _d_Chr_Cccc.Atar + ( Chr_Ccc.Ata 
               * Chr_Ccc.Ttr + Chr_Ccc.Ara * Chr_Ccc.Rtr + Chr_Ccc.Aaa 
               * Chr_Ccc.Atr + Chr_Ccc.Apa * Chr_Ccc.Ptr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tta + Chr_Ccc.Arr * Chr_Ccc.Rta + Chr_Ccc.Aar 
               * Chr_Ccc.Ata + Chr_Ccc.Apr * Chr_Ccc.Pta );
 R_Cccc.Ataa = _d_Chr_Cccc.Ataa - _d_Chr_Cccc.Ataa + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tta + Chr_Ccc.Ara * Chr_Ccc.Rta + Chr_Ccc.Aaa 
               * Chr_Ccc.Ata + Chr_Ccc.Apa * Chr_Ccc.Pta ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tta + Chr_Ccc.Ara * Chr_Ccc.Rta + Chr_Ccc.Aaa 
               * Chr_Ccc.Ata + Chr_Ccc.Apa * Chr_Ccc.Pta );
 R_Cccc.Atap = _d_Chr_Cccc.Atpa - _d_Chr_Cccc.Atap + ( Chr_Ccc.Ata 
               * Chr_Ccc.Ttp + Chr_Ccc.Ara * Chr_Ccc.Rtp + Chr_Ccc.Aaa 
               * Chr_Ccc.Atp + Chr_Ccc.Apa * Chr_Ccc.Ptp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tta + Chr_Ccc.Arp * Chr_Ccc.Rta + Chr_Ccc.Aap 
               * Chr_Ccc.Ata + Chr_Ccc.App * Chr_Ccc.Pta );
 R_Cccc.Atpt = _d_Chr_Cccc.Attp - _d_Chr_Cccc.Atpt + ( Chr_Ccc.Atp 
               * Chr_Ccc.Ttt + Chr_Ccc.Arp * Chr_Ccc.Rtt + Chr_Ccc.Aap 
               * Chr_Ccc.Att + Chr_Ccc.App * Chr_Ccc.Ptt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Ttp + Chr_Ccc.Art * Chr_Ccc.Rtp + Chr_Ccc.Aat 
               * Chr_Ccc.Atp + Chr_Ccc.Apt * Chr_Ccc.Ptp );
 R_Cccc.Atpr = _d_Chr_Cccc.Atrp - _d_Chr_Cccc.Atpr + ( Chr_Ccc.Atp 
               * Chr_Ccc.Ttr + Chr_Ccc.Arp * Chr_Ccc.Rtr + Chr_Ccc.Aap 
               * Chr_Ccc.Atr + Chr_Ccc.App * Chr_Ccc.Ptr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Ttp + Chr_Ccc.Arr * Chr_Ccc.Rtp + Chr_Ccc.Aar 
               * Chr_Ccc.Atp + Chr_Ccc.Apr * Chr_Ccc.Ptp );
 R_Cccc.Atpa = _d_Chr_Cccc.Atap - _d_Chr_Cccc.Atpa + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tta + Chr_Ccc.Arp * Chr_Ccc.Rta + Chr_Ccc.Aap 
               * Chr_Ccc.Ata + Chr_Ccc.App * Chr_Ccc.Pta ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Ttp + Chr_Ccc.Ara * Chr_Ccc.Rtp + Chr_Ccc.Aaa 
               * Chr_Ccc.Atp + Chr_Ccc.Apa * Chr_Ccc.Ptp );
 R_Cccc.Atpp = _d_Chr_Cccc.Atpp - _d_Chr_Cccc.Atpp + ( Chr_Ccc.Atp 
               * Chr_Ccc.Ttp + Chr_Ccc.Arp * Chr_Ccc.Rtp + Chr_Ccc.Aap 
               * Chr_Ccc.Atp + Chr_Ccc.App * Chr_Ccc.Ptp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Ttp + Chr_Ccc.Arp * Chr_Ccc.Rtp + Chr_Ccc.Aap 
               * Chr_Ccc.Atp + Chr_Ccc.App * Chr_Ccc.Ptp );
 R_Cccc.Artt = _d_Chr_Cccc.Artt - _d_Chr_Cccc.Artt + ( Chr_Ccc.Att 
               * Chr_Ccc.Trt + Chr_Ccc.Art * Chr_Ccc.Rrt + Chr_Ccc.Aat 
               * Chr_Ccc.Art + Chr_Ccc.Apt * Chr_Ccc.Prt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Trt + Chr_Ccc.Art * Chr_Ccc.Rrt + Chr_Ccc.Aat 
               * Chr_Ccc.Art + Chr_Ccc.Apt * Chr_Ccc.Prt );
 R_Cccc.Artr = _d_Chr_Cccc.Arrt - _d_Chr_Cccc.Artr + ( Chr_Ccc.Att 
               * Chr_Ccc.Trr + Chr_Ccc.Art * Chr_Ccc.Rrr + Chr_Ccc.Aat 
               * Chr_Ccc.Arr + Chr_Ccc.Apt * Chr_Ccc.Prr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Trt + Chr_Ccc.Arr * Chr_Ccc.Rrt + Chr_Ccc.Aar 
               * Chr_Ccc.Art + Chr_Ccc.Apr * Chr_Ccc.Prt );
 R_Cccc.Arta = _d_Chr_Cccc.Arat - _d_Chr_Cccc.Arta + ( Chr_Ccc.Att 
               * Chr_Ccc.Tra + Chr_Ccc.Art * Chr_Ccc.Rra + Chr_Ccc.Aat 
               * Chr_Ccc.Ara + Chr_Ccc.Apt * Chr_Ccc.Pra ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Trt + Chr_Ccc.Ara * Chr_Ccc.Rrt + Chr_Ccc.Aaa 
               * Chr_Ccc.Art + Chr_Ccc.Apa * Chr_Ccc.Prt );
 R_Cccc.Artp = _d_Chr_Cccc.Arpt - _d_Chr_Cccc.Artp + ( Chr_Ccc.Att 
               * Chr_Ccc.Trp + Chr_Ccc.Art * Chr_Ccc.Rrp + Chr_Ccc.Aat 
               * Chr_Ccc.Arp + Chr_Ccc.Apt * Chr_Ccc.Prp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Trt + Chr_Ccc.Arp * Chr_Ccc.Rrt + Chr_Ccc.Aap 
               * Chr_Ccc.Art + Chr_Ccc.App * Chr_Ccc.Prt );
 R_Cccc.Arrt = _d_Chr_Cccc.Artr - _d_Chr_Cccc.Arrt + ( Chr_Ccc.Atr 
               * Chr_Ccc.Trt + Chr_Ccc.Arr * Chr_Ccc.Rrt + Chr_Ccc.Aar 
               * Chr_Ccc.Art + Chr_Ccc.Apr * Chr_Ccc.Prt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Trr + Chr_Ccc.Art * Chr_Ccc.Rrr + Chr_Ccc.Aat 
               * Chr_Ccc.Arr + Chr_Ccc.Apt * Chr_Ccc.Prr );
 R_Cccc.Arrr = _d_Chr_Cccc.Arrr - _d_Chr_Cccc.Arrr + ( Chr_Ccc.Atr 
               * Chr_Ccc.Trr + Chr_Ccc.Arr * Chr_Ccc.Rrr + Chr_Ccc.Aar 
               * Chr_Ccc.Arr + Chr_Ccc.Apr * Chr_Ccc.Prr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Trr + Chr_Ccc.Arr * Chr_Ccc.Rrr + Chr_Ccc.Aar 
               * Chr_Ccc.Arr + Chr_Ccc.Apr * Chr_Ccc.Prr );
 R_Cccc.Arra = _d_Chr_Cccc.Arar - _d_Chr_Cccc.Arra + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tra + Chr_Ccc.Arr * Chr_Ccc.Rra + Chr_Ccc.Aar 
               * Chr_Ccc.Ara + Chr_Ccc.Apr * Chr_Ccc.Pra ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Trr + Chr_Ccc.Ara * Chr_Ccc.Rrr + Chr_Ccc.Aaa 
               * Chr_Ccc.Arr + Chr_Ccc.Apa * Chr_Ccc.Prr );
 R_Cccc.Arrp = _d_Chr_Cccc.Arpr - _d_Chr_Cccc.Arrp + ( Chr_Ccc.Atr 
               * Chr_Ccc.Trp + Chr_Ccc.Arr * Chr_Ccc.Rrp + Chr_Ccc.Aar 
               * Chr_Ccc.Arp + Chr_Ccc.Apr * Chr_Ccc.Prp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Trr + Chr_Ccc.Arp * Chr_Ccc.Rrr + Chr_Ccc.Aap 
               * Chr_Ccc.Arr + Chr_Ccc.App * Chr_Ccc.Prr );
 R_Cccc.Arat = _d_Chr_Cccc.Arta - _d_Chr_Cccc.Arat + ( Chr_Ccc.Ata 
               * Chr_Ccc.Trt + Chr_Ccc.Ara * Chr_Ccc.Rrt + Chr_Ccc.Aaa 
               * Chr_Ccc.Art + Chr_Ccc.Apa * Chr_Ccc.Prt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tra + Chr_Ccc.Art * Chr_Ccc.Rra + Chr_Ccc.Aat 
               * Chr_Ccc.Ara + Chr_Ccc.Apt * Chr_Ccc.Pra );
 R_Cccc.Arar = _d_Chr_Cccc.Arra - _d_Chr_Cccc.Arar + ( Chr_Ccc.Ata 
               * Chr_Ccc.Trr + Chr_Ccc.Ara * Chr_Ccc.Rrr + Chr_Ccc.Aaa 
               * Chr_Ccc.Arr + Chr_Ccc.Apa * Chr_Ccc.Prr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tra + Chr_Ccc.Arr * Chr_Ccc.Rra + Chr_Ccc.Aar 
               * Chr_Ccc.Ara + Chr_Ccc.Apr * Chr_Ccc.Pra );
 R_Cccc.Araa = _d_Chr_Cccc.Araa - _d_Chr_Cccc.Araa + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tra + Chr_Ccc.Ara * Chr_Ccc.Rra + Chr_Ccc.Aaa 
               * Chr_Ccc.Ara + Chr_Ccc.Apa * Chr_Ccc.Pra ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tra + Chr_Ccc.Ara * Chr_Ccc.Rra + Chr_Ccc.Aaa 
               * Chr_Ccc.Ara + Chr_Ccc.Apa * Chr_Ccc.Pra );
 R_Cccc.Arap = _d_Chr_Cccc.Arpa - _d_Chr_Cccc.Arap + ( Chr_Ccc.Ata 
               * Chr_Ccc.Trp + Chr_Ccc.Ara * Chr_Ccc.Rrp + Chr_Ccc.Aaa 
               * Chr_Ccc.Arp + Chr_Ccc.Apa * Chr_Ccc.Prp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tra + Chr_Ccc.Arp * Chr_Ccc.Rra + Chr_Ccc.Aap 
               * Chr_Ccc.Ara + Chr_Ccc.App * Chr_Ccc.Pra );
 R_Cccc.Arpt = _d_Chr_Cccc.Artp - _d_Chr_Cccc.Arpt + ( Chr_Ccc.Atp 
               * Chr_Ccc.Trt + Chr_Ccc.Arp * Chr_Ccc.Rrt + Chr_Ccc.Aap 
               * Chr_Ccc.Art + Chr_Ccc.App * Chr_Ccc.Prt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Trp + Chr_Ccc.Art * Chr_Ccc.Rrp + Chr_Ccc.Aat 
               * Chr_Ccc.Arp + Chr_Ccc.Apt * Chr_Ccc.Prp );
 R_Cccc.Arpr = _d_Chr_Cccc.Arrp - _d_Chr_Cccc.Arpr + ( Chr_Ccc.Atp 
               * Chr_Ccc.Trr + Chr_Ccc.Arp * Chr_Ccc.Rrr + Chr_Ccc.Aap 
               * Chr_Ccc.Arr + Chr_Ccc.App * Chr_Ccc.Prr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Trp + Chr_Ccc.Arr * Chr_Ccc.Rrp + Chr_Ccc.Aar 
               * Chr_Ccc.Arp + Chr_Ccc.Apr * Chr_Ccc.Prp );
 R_Cccc.Arpa = _d_Chr_Cccc.Arap - _d_Chr_Cccc.Arpa + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tra + Chr_Ccc.Arp * Chr_Ccc.Rra + Chr_Ccc.Aap 
               * Chr_Ccc.Ara + Chr_Ccc.App * Chr_Ccc.Pra ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Trp + Chr_Ccc.Ara * Chr_Ccc.Rrp + Chr_Ccc.Aaa 
               * Chr_Ccc.Arp + Chr_Ccc.Apa * Chr_Ccc.Prp );
 R_Cccc.Arpp = _d_Chr_Cccc.Arpp - _d_Chr_Cccc.Arpp + ( Chr_Ccc.Atp 
               * Chr_Ccc.Trp + Chr_Ccc.Arp * Chr_Ccc.Rrp + Chr_Ccc.Aap 
               * Chr_Ccc.Arp + Chr_Ccc.App * Chr_Ccc.Prp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Trp + Chr_Ccc.Arp * Chr_Ccc.Rrp + Chr_Ccc.Aap 
               * Chr_Ccc.Arp + Chr_Ccc.App * Chr_Ccc.Prp );
 R_Cccc.Aatt = _d_Chr_Cccc.Aatt - _d_Chr_Cccc.Aatt + ( Chr_Ccc.Att 
               * Chr_Ccc.Tat + Chr_Ccc.Art * Chr_Ccc.Rat + Chr_Ccc.Aat 
               * Chr_Ccc.Aat + Chr_Ccc.Apt * Chr_Ccc.Pat ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tat + Chr_Ccc.Art * Chr_Ccc.Rat + Chr_Ccc.Aat 
               * Chr_Ccc.Aat + Chr_Ccc.Apt * Chr_Ccc.Pat );
 R_Cccc.Aatr = _d_Chr_Cccc.Aart - _d_Chr_Cccc.Aatr + ( Chr_Ccc.Att 
               * Chr_Ccc.Tar + Chr_Ccc.Art * Chr_Ccc.Rar + Chr_Ccc.Aat 
               * Chr_Ccc.Aar + Chr_Ccc.Apt * Chr_Ccc.Par ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tat + Chr_Ccc.Arr * Chr_Ccc.Rat + Chr_Ccc.Aar 
               * Chr_Ccc.Aat + Chr_Ccc.Apr * Chr_Ccc.Pat );
 R_Cccc.Aata = _d_Chr_Cccc.Aaat - _d_Chr_Cccc.Aata + ( Chr_Ccc.Att 
               * Chr_Ccc.Taa + Chr_Ccc.Art * Chr_Ccc.Raa + Chr_Ccc.Aat 
               * Chr_Ccc.Aaa + Chr_Ccc.Apt * Chr_Ccc.Paa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tat + Chr_Ccc.Ara * Chr_Ccc.Rat + Chr_Ccc.Aaa 
               * Chr_Ccc.Aat + Chr_Ccc.Apa * Chr_Ccc.Pat );
 R_Cccc.Aatp = _d_Chr_Cccc.Aapt - _d_Chr_Cccc.Aatp + ( Chr_Ccc.Att 
               * Chr_Ccc.Tap + Chr_Ccc.Art * Chr_Ccc.Rap + Chr_Ccc.Aat 
               * Chr_Ccc.Aap + Chr_Ccc.Apt * Chr_Ccc.Pap ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tat + Chr_Ccc.Arp * Chr_Ccc.Rat + Chr_Ccc.Aap 
               * Chr_Ccc.Aat + Chr_Ccc.App * Chr_Ccc.Pat );
 R_Cccc.Aart = _d_Chr_Cccc.Aatr - _d_Chr_Cccc.Aart + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tat + Chr_Ccc.Arr * Chr_Ccc.Rat + Chr_Ccc.Aar 
               * Chr_Ccc.Aat + Chr_Ccc.Apr * Chr_Ccc.Pat ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tar + Chr_Ccc.Art * Chr_Ccc.Rar + Chr_Ccc.Aat 
               * Chr_Ccc.Aar + Chr_Ccc.Apt * Chr_Ccc.Par );
 R_Cccc.Aarr = _d_Chr_Cccc.Aarr - _d_Chr_Cccc.Aarr + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tar + Chr_Ccc.Arr * Chr_Ccc.Rar + Chr_Ccc.Aar 
               * Chr_Ccc.Aar + Chr_Ccc.Apr * Chr_Ccc.Par ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tar + Chr_Ccc.Arr * Chr_Ccc.Rar + Chr_Ccc.Aar 
               * Chr_Ccc.Aar + Chr_Ccc.Apr * Chr_Ccc.Par );
 R_Cccc.Aara = _d_Chr_Cccc.Aaar - _d_Chr_Cccc.Aara + ( Chr_Ccc.Atr 
               * Chr_Ccc.Taa + Chr_Ccc.Arr * Chr_Ccc.Raa + Chr_Ccc.Aar 
               * Chr_Ccc.Aaa + Chr_Ccc.Apr * Chr_Ccc.Paa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tar + Chr_Ccc.Ara * Chr_Ccc.Rar + Chr_Ccc.Aaa 
               * Chr_Ccc.Aar + Chr_Ccc.Apa * Chr_Ccc.Par );
 R_Cccc.Aarp = _d_Chr_Cccc.Aapr - _d_Chr_Cccc.Aarp + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tap + Chr_Ccc.Arr * Chr_Ccc.Rap + Chr_Ccc.Aar 
               * Chr_Ccc.Aap + Chr_Ccc.Apr * Chr_Ccc.Pap ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tar + Chr_Ccc.Arp * Chr_Ccc.Rar + Chr_Ccc.Aap 
               * Chr_Ccc.Aar + Chr_Ccc.App * Chr_Ccc.Par );
 R_Cccc.Aaat = _d_Chr_Cccc.Aata - _d_Chr_Cccc.Aaat + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tat + Chr_Ccc.Ara * Chr_Ccc.Rat + Chr_Ccc.Aaa 
               * Chr_Ccc.Aat + Chr_Ccc.Apa * Chr_Ccc.Pat ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Taa + Chr_Ccc.Art * Chr_Ccc.Raa + Chr_Ccc.Aat 
               * Chr_Ccc.Aaa + Chr_Ccc.Apt * Chr_Ccc.Paa );
 R_Cccc.Aaar = _d_Chr_Cccc.Aara - _d_Chr_Cccc.Aaar + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tar + Chr_Ccc.Ara * Chr_Ccc.Rar + Chr_Ccc.Aaa 
               * Chr_Ccc.Aar + Chr_Ccc.Apa * Chr_Ccc.Par ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Taa + Chr_Ccc.Arr * Chr_Ccc.Raa + Chr_Ccc.Aar 
               * Chr_Ccc.Aaa + Chr_Ccc.Apr * Chr_Ccc.Paa );
 R_Cccc.Aaaa = _d_Chr_Cccc.Aaaa - _d_Chr_Cccc.Aaaa + ( Chr_Ccc.Ata 
               * Chr_Ccc.Taa + Chr_Ccc.Ara * Chr_Ccc.Raa + Chr_Ccc.Aaa 
               * Chr_Ccc.Aaa + Chr_Ccc.Apa * Chr_Ccc.Paa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Taa + Chr_Ccc.Ara * Chr_Ccc.Raa + Chr_Ccc.Aaa 
               * Chr_Ccc.Aaa + Chr_Ccc.Apa * Chr_Ccc.Paa );
 R_Cccc.Aaap = _d_Chr_Cccc.Aapa - _d_Chr_Cccc.Aaap + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tap + Chr_Ccc.Ara * Chr_Ccc.Rap + Chr_Ccc.Aaa 
               * Chr_Ccc.Aap + Chr_Ccc.Apa * Chr_Ccc.Pap ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Taa + Chr_Ccc.Arp * Chr_Ccc.Raa + Chr_Ccc.Aap 
               * Chr_Ccc.Aaa + Chr_Ccc.App * Chr_Ccc.Paa );
 R_Cccc.Aapt = _d_Chr_Cccc.Aatp - _d_Chr_Cccc.Aapt + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tat + Chr_Ccc.Arp * Chr_Ccc.Rat + Chr_Ccc.Aap 
               * Chr_Ccc.Aat + Chr_Ccc.App * Chr_Ccc.Pat ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tap + Chr_Ccc.Art * Chr_Ccc.Rap + Chr_Ccc.Aat 
               * Chr_Ccc.Aap + Chr_Ccc.Apt * Chr_Ccc.Pap );
 R_Cccc.Aapr = _d_Chr_Cccc.Aarp - _d_Chr_Cccc.Aapr + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tar + Chr_Ccc.Arp * Chr_Ccc.Rar + Chr_Ccc.Aap 
               * Chr_Ccc.Aar + Chr_Ccc.App * Chr_Ccc.Par ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tap + Chr_Ccc.Arr * Chr_Ccc.Rap + Chr_Ccc.Aar 
               * Chr_Ccc.Aap + Chr_Ccc.Apr * Chr_Ccc.Pap );
 R_Cccc.Aapa = _d_Chr_Cccc.Aaap - _d_Chr_Cccc.Aapa + ( Chr_Ccc.Atp 
               * Chr_Ccc.Taa + Chr_Ccc.Arp * Chr_Ccc.Raa + Chr_Ccc.Aap 
               * Chr_Ccc.Aaa + Chr_Ccc.App * Chr_Ccc.Paa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tap + Chr_Ccc.Ara * Chr_Ccc.Rap + Chr_Ccc.Aaa 
               * Chr_Ccc.Aap + Chr_Ccc.Apa * Chr_Ccc.Pap );
 R_Cccc.Aapp = _d_Chr_Cccc.Aapp - _d_Chr_Cccc.Aapp + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tap + Chr_Ccc.Arp * Chr_Ccc.Rap + Chr_Ccc.Aap 
               * Chr_Ccc.Aap + Chr_Ccc.App * Chr_Ccc.Pap ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tap + Chr_Ccc.Arp * Chr_Ccc.Rap + Chr_Ccc.Aap 
               * Chr_Ccc.Aap + Chr_Ccc.App * Chr_Ccc.Pap );
 R_Cccc.Aptt = _d_Chr_Cccc.Aptt - _d_Chr_Cccc.Aptt + ( Chr_Ccc.Att 
               * Chr_Ccc.Tpt + Chr_Ccc.Art * Chr_Ccc.Rpt + Chr_Ccc.Aat 
               * Chr_Ccc.Apt + Chr_Ccc.Apt * Chr_Ccc.Ppt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tpt + Chr_Ccc.Art * Chr_Ccc.Rpt + Chr_Ccc.Aat 
               * Chr_Ccc.Apt + Chr_Ccc.Apt * Chr_Ccc.Ppt );
 R_Cccc.Aptr = _d_Chr_Cccc.Aprt - _d_Chr_Cccc.Aptr + ( Chr_Ccc.Att 
               * Chr_Ccc.Tpr + Chr_Ccc.Art * Chr_Ccc.Rpr + Chr_Ccc.Aat 
               * Chr_Ccc.Apr + Chr_Ccc.Apt * Chr_Ccc.Ppr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpt + Chr_Ccc.Arr * Chr_Ccc.Rpt + Chr_Ccc.Aar 
               * Chr_Ccc.Apt + Chr_Ccc.Apr * Chr_Ccc.Ppt );
 R_Cccc.Apta = _d_Chr_Cccc.Apat - _d_Chr_Cccc.Apta + ( Chr_Ccc.Att 
               * Chr_Ccc.Tpa + Chr_Ccc.Art * Chr_Ccc.Rpa + Chr_Ccc.Aat 
               * Chr_Ccc.Apa + Chr_Ccc.Apt * Chr_Ccc.Ppa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpt + Chr_Ccc.Ara * Chr_Ccc.Rpt + Chr_Ccc.Aaa 
               * Chr_Ccc.Apt + Chr_Ccc.Apa * Chr_Ccc.Ppt );
 R_Cccc.Aptp = _d_Chr_Cccc.Appt - _d_Chr_Cccc.Aptp + ( Chr_Ccc.Att 
               * Chr_Ccc.Tpp + Chr_Ccc.Art * Chr_Ccc.Rpp + Chr_Ccc.Aat 
               * Chr_Ccc.App + Chr_Ccc.Apt * Chr_Ccc.Ppp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpt + Chr_Ccc.Arp * Chr_Ccc.Rpt + Chr_Ccc.Aap 
               * Chr_Ccc.Apt + Chr_Ccc.App * Chr_Ccc.Ppt );
 R_Cccc.Aprt = _d_Chr_Cccc.Aptr - _d_Chr_Cccc.Aprt + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpt + Chr_Ccc.Arr * Chr_Ccc.Rpt + Chr_Ccc.Aar 
               * Chr_Ccc.Apt + Chr_Ccc.Apr * Chr_Ccc.Ppt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tpr + Chr_Ccc.Art * Chr_Ccc.Rpr + Chr_Ccc.Aat 
               * Chr_Ccc.Apr + Chr_Ccc.Apt * Chr_Ccc.Ppr );
 R_Cccc.Aprr = _d_Chr_Cccc.Aprr - _d_Chr_Cccc.Aprr + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpr + Chr_Ccc.Arr * Chr_Ccc.Rpr + Chr_Ccc.Aar 
               * Chr_Ccc.Apr + Chr_Ccc.Apr * Chr_Ccc.Ppr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpr + Chr_Ccc.Arr * Chr_Ccc.Rpr + Chr_Ccc.Aar 
               * Chr_Ccc.Apr + Chr_Ccc.Apr * Chr_Ccc.Ppr );
 R_Cccc.Apra = _d_Chr_Cccc.Apar - _d_Chr_Cccc.Apra + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpa + Chr_Ccc.Arr * Chr_Ccc.Rpa + Chr_Ccc.Aar 
               * Chr_Ccc.Apa + Chr_Ccc.Apr * Chr_Ccc.Ppa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpr + Chr_Ccc.Ara * Chr_Ccc.Rpr + Chr_Ccc.Aaa 
               * Chr_Ccc.Apr + Chr_Ccc.Apa * Chr_Ccc.Ppr );
 R_Cccc.Aprp = _d_Chr_Cccc.Appr - _d_Chr_Cccc.Aprp + ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpp + Chr_Ccc.Arr * Chr_Ccc.Rpp + Chr_Ccc.Aar 
               * Chr_Ccc.App + Chr_Ccc.Apr * Chr_Ccc.Ppp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpr + Chr_Ccc.Arp * Chr_Ccc.Rpr + Chr_Ccc.Aap 
               * Chr_Ccc.Apr + Chr_Ccc.App * Chr_Ccc.Ppr );
 R_Cccc.Apat = _d_Chr_Cccc.Apta - _d_Chr_Cccc.Apat + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpt + Chr_Ccc.Ara * Chr_Ccc.Rpt + Chr_Ccc.Aaa 
               * Chr_Ccc.Apt + Chr_Ccc.Apa * Chr_Ccc.Ppt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tpa + Chr_Ccc.Art * Chr_Ccc.Rpa + Chr_Ccc.Aat 
               * Chr_Ccc.Apa + Chr_Ccc.Apt * Chr_Ccc.Ppa );
 R_Cccc.Apar = _d_Chr_Cccc.Apra - _d_Chr_Cccc.Apar + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpr + Chr_Ccc.Ara * Chr_Ccc.Rpr + Chr_Ccc.Aaa 
               * Chr_Ccc.Apr + Chr_Ccc.Apa * Chr_Ccc.Ppr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpa + Chr_Ccc.Arr * Chr_Ccc.Rpa + Chr_Ccc.Aar 
               * Chr_Ccc.Apa + Chr_Ccc.Apr * Chr_Ccc.Ppa );
 R_Cccc.Apaa = _d_Chr_Cccc.Apaa - _d_Chr_Cccc.Apaa + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpa + Chr_Ccc.Ara * Chr_Ccc.Rpa + Chr_Ccc.Aaa 
               * Chr_Ccc.Apa + Chr_Ccc.Apa * Chr_Ccc.Ppa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpa + Chr_Ccc.Ara * Chr_Ccc.Rpa + Chr_Ccc.Aaa 
               * Chr_Ccc.Apa + Chr_Ccc.Apa * Chr_Ccc.Ppa );
 R_Cccc.Apap = _d_Chr_Cccc.Appa - _d_Chr_Cccc.Apap + ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpp + Chr_Ccc.Ara * Chr_Ccc.Rpp + Chr_Ccc.Aaa 
               * Chr_Ccc.App + Chr_Ccc.Apa * Chr_Ccc.Ppp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpa + Chr_Ccc.Arp * Chr_Ccc.Rpa + Chr_Ccc.Aap 
               * Chr_Ccc.Apa + Chr_Ccc.App * Chr_Ccc.Ppa );
 R_Cccc.Appt = _d_Chr_Cccc.Aptp - _d_Chr_Cccc.Appt + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpt + Chr_Ccc.Arp * Chr_Ccc.Rpt + Chr_Ccc.Aap 
               * Chr_Ccc.Apt + Chr_Ccc.App * Chr_Ccc.Ppt ) - ( Chr_Ccc.Att 
               * Chr_Ccc.Tpp + Chr_Ccc.Art * Chr_Ccc.Rpp + Chr_Ccc.Aat 
               * Chr_Ccc.App + Chr_Ccc.Apt * Chr_Ccc.Ppp );
 R_Cccc.Appr = _d_Chr_Cccc.Aprp - _d_Chr_Cccc.Appr + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpr + Chr_Ccc.Arp * Chr_Ccc.Rpr + Chr_Ccc.Aap 
               * Chr_Ccc.Apr + Chr_Ccc.App * Chr_Ccc.Ppr ) - ( Chr_Ccc.Atr 
               * Chr_Ccc.Tpp + Chr_Ccc.Arr * Chr_Ccc.Rpp + Chr_Ccc.Aar 
               * Chr_Ccc.App + Chr_Ccc.Apr * Chr_Ccc.Ppp );
 R_Cccc.Appa = _d_Chr_Cccc.Apap - _d_Chr_Cccc.Appa + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpa + Chr_Ccc.Arp * Chr_Ccc.Rpa + Chr_Ccc.Aap 
               * Chr_Ccc.Apa + Chr_Ccc.App * Chr_Ccc.Ppa ) - ( Chr_Ccc.Ata 
               * Chr_Ccc.Tpp + Chr_Ccc.Ara * Chr_Ccc.Rpp + Chr_Ccc.Aaa 
               * Chr_Ccc.App + Chr_Ccc.Apa * Chr_Ccc.Ppp );
 R_Cccc.Appp = _d_Chr_Cccc.Appp - _d_Chr_Cccc.Appp + ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpp + Chr_Ccc.Arp * Chr_Ccc.Rpp + Chr_Ccc.Aap 
               * Chr_Ccc.App + Chr_Ccc.App * Chr_Ccc.Ppp ) - ( Chr_Ccc.Atp 
               * Chr_Ccc.Tpp + Chr_Ccc.Arp * Chr_Ccc.Rpp + Chr_Ccc.Aap 
               * Chr_Ccc.App + Chr_Ccc.App * Chr_Ccc.Ppp );
 R_Cccc.Pttt = _d_Chr_Cccc.Pttt - _d_Chr_Cccc.Pttt + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Ttt + Chr_Ccc.Prt * Chr_Ccc.Rtt + Chr_Ccc.Pat 
               * Chr_Ccc.Att + Chr_Ccc.Ppt * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Ttt + Chr_Ccc.Prt * Chr_Ccc.Rtt + Chr_Ccc.Pat 
               * Chr_Ccc.Att + Chr_Ccc.Ppt * Chr_Ccc.Ptt );
 R_Cccc.Pttr = _d_Chr_Cccc.Ptrt - _d_Chr_Cccc.Pttr + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Ttr + Chr_Ccc.Prt * Chr_Ccc.Rtr + Chr_Ccc.Pat 
               * Chr_Ccc.Atr + Chr_Ccc.Ppt * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Ttt + Chr_Ccc.Prr * Chr_Ccc.Rtt + Chr_Ccc.Par 
               * Chr_Ccc.Att + Chr_Ccc.Ppr * Chr_Ccc.Ptt );
 R_Cccc.Ptta = _d_Chr_Cccc.Ptat - _d_Chr_Cccc.Ptta + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tta + Chr_Ccc.Prt * Chr_Ccc.Rta + Chr_Ccc.Pat 
               * Chr_Ccc.Ata + Chr_Ccc.Ppt * Chr_Ccc.Pta ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Ttt + Chr_Ccc.Pra * Chr_Ccc.Rtt + Chr_Ccc.Paa 
               * Chr_Ccc.Att + Chr_Ccc.Ppa * Chr_Ccc.Ptt );
 R_Cccc.Pttp = _d_Chr_Cccc.Ptpt - _d_Chr_Cccc.Pttp + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Ttp + Chr_Ccc.Prt * Chr_Ccc.Rtp + Chr_Ccc.Pat 
               * Chr_Ccc.Atp + Chr_Ccc.Ppt * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Ttt + Chr_Ccc.Prp * Chr_Ccc.Rtt + Chr_Ccc.Pap 
               * Chr_Ccc.Att + Chr_Ccc.Ppp * Chr_Ccc.Ptt );
 R_Cccc.Ptrt = _d_Chr_Cccc.Pttr - _d_Chr_Cccc.Ptrt + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Ttt + Chr_Ccc.Prr * Chr_Ccc.Rtt + Chr_Ccc.Par 
               * Chr_Ccc.Att + Chr_Ccc.Ppr * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Ttr + Chr_Ccc.Prt * Chr_Ccc.Rtr + Chr_Ccc.Pat 
               * Chr_Ccc.Atr + Chr_Ccc.Ppt * Chr_Ccc.Ptr );
 R_Cccc.Ptrr = _d_Chr_Cccc.Ptrr - _d_Chr_Cccc.Ptrr + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Ttr + Chr_Ccc.Prr * Chr_Ccc.Rtr + Chr_Ccc.Par 
               * Chr_Ccc.Atr + Chr_Ccc.Ppr * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Ttr + Chr_Ccc.Prr * Chr_Ccc.Rtr + Chr_Ccc.Par 
               * Chr_Ccc.Atr + Chr_Ccc.Ppr * Chr_Ccc.Ptr );
 R_Cccc.Ptra = _d_Chr_Cccc.Ptar - _d_Chr_Cccc.Ptra + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tta + Chr_Ccc.Prr * Chr_Ccc.Rta + Chr_Ccc.Par 
               * Chr_Ccc.Ata + Chr_Ccc.Ppr * Chr_Ccc.Pta ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Ttr + Chr_Ccc.Pra * Chr_Ccc.Rtr + Chr_Ccc.Paa 
               * Chr_Ccc.Atr + Chr_Ccc.Ppa * Chr_Ccc.Ptr );
 R_Cccc.Ptrp = _d_Chr_Cccc.Ptpr - _d_Chr_Cccc.Ptrp + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Ttp + Chr_Ccc.Prr * Chr_Ccc.Rtp + Chr_Ccc.Par 
               * Chr_Ccc.Atp + Chr_Ccc.Ppr * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Ttr + Chr_Ccc.Prp * Chr_Ccc.Rtr + Chr_Ccc.Pap 
               * Chr_Ccc.Atr + Chr_Ccc.Ppp * Chr_Ccc.Ptr );
 R_Cccc.Ptat = _d_Chr_Cccc.Ptta - _d_Chr_Cccc.Ptat + ( Chr_Ccc.Pta 
               * Chr_Ccc.Ttt + Chr_Ccc.Pra * Chr_Ccc.Rtt + Chr_Ccc.Paa 
               * Chr_Ccc.Att + Chr_Ccc.Ppa * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tta + Chr_Ccc.Prt * Chr_Ccc.Rta + Chr_Ccc.Pat 
               * Chr_Ccc.Ata + Chr_Ccc.Ppt * Chr_Ccc.Pta );
 R_Cccc.Ptar = _d_Chr_Cccc.Ptra - _d_Chr_Cccc.Ptar + ( Chr_Ccc.Pta 
               * Chr_Ccc.Ttr + Chr_Ccc.Pra * Chr_Ccc.Rtr + Chr_Ccc.Paa 
               * Chr_Ccc.Atr + Chr_Ccc.Ppa * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tta + Chr_Ccc.Prr * Chr_Ccc.Rta + Chr_Ccc.Par 
               * Chr_Ccc.Ata + Chr_Ccc.Ppr * Chr_Ccc.Pta );
 R_Cccc.Ptaa = _d_Chr_Cccc.Ptaa - _d_Chr_Cccc.Ptaa + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tta + Chr_Ccc.Pra * Chr_Ccc.Rta + Chr_Ccc.Paa 
               * Chr_Ccc.Ata + Chr_Ccc.Ppa * Chr_Ccc.Pta ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tta + Chr_Ccc.Pra * Chr_Ccc.Rta + Chr_Ccc.Paa 
               * Chr_Ccc.Ata + Chr_Ccc.Ppa * Chr_Ccc.Pta );
 R_Cccc.Ptap = _d_Chr_Cccc.Ptpa - _d_Chr_Cccc.Ptap + ( Chr_Ccc.Pta 
               * Chr_Ccc.Ttp + Chr_Ccc.Pra * Chr_Ccc.Rtp + Chr_Ccc.Paa 
               * Chr_Ccc.Atp + Chr_Ccc.Ppa * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tta + Chr_Ccc.Prp * Chr_Ccc.Rta + Chr_Ccc.Pap 
               * Chr_Ccc.Ata + Chr_Ccc.Ppp * Chr_Ccc.Pta );
 R_Cccc.Ptpt = _d_Chr_Cccc.Pttp - _d_Chr_Cccc.Ptpt + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Ttt + Chr_Ccc.Prp * Chr_Ccc.Rtt + Chr_Ccc.Pap 
               * Chr_Ccc.Att + Chr_Ccc.Ppp * Chr_Ccc.Ptt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Ttp + Chr_Ccc.Prt * Chr_Ccc.Rtp + Chr_Ccc.Pat 
               * Chr_Ccc.Atp + Chr_Ccc.Ppt * Chr_Ccc.Ptp );
 R_Cccc.Ptpr = _d_Chr_Cccc.Ptrp - _d_Chr_Cccc.Ptpr + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Ttr + Chr_Ccc.Prp * Chr_Ccc.Rtr + Chr_Ccc.Pap 
               * Chr_Ccc.Atr + Chr_Ccc.Ppp * Chr_Ccc.Ptr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Ttp + Chr_Ccc.Prr * Chr_Ccc.Rtp + Chr_Ccc.Par 
               * Chr_Ccc.Atp + Chr_Ccc.Ppr * Chr_Ccc.Ptp );
 R_Cccc.Ptpa = _d_Chr_Cccc.Ptap - _d_Chr_Cccc.Ptpa + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tta + Chr_Ccc.Prp * Chr_Ccc.Rta + Chr_Ccc.Pap 
               * Chr_Ccc.Ata + Chr_Ccc.Ppp * Chr_Ccc.Pta ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Ttp + Chr_Ccc.Pra * Chr_Ccc.Rtp + Chr_Ccc.Paa 
               * Chr_Ccc.Atp + Chr_Ccc.Ppa * Chr_Ccc.Ptp );
 R_Cccc.Ptpp = _d_Chr_Cccc.Ptpp - _d_Chr_Cccc.Ptpp + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Ttp + Chr_Ccc.Prp * Chr_Ccc.Rtp + Chr_Ccc.Pap 
               * Chr_Ccc.Atp + Chr_Ccc.Ppp * Chr_Ccc.Ptp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Ttp + Chr_Ccc.Prp * Chr_Ccc.Rtp + Chr_Ccc.Pap 
               * Chr_Ccc.Atp + Chr_Ccc.Ppp * Chr_Ccc.Ptp );
 R_Cccc.Prtt = _d_Chr_Cccc.Prtt - _d_Chr_Cccc.Prtt + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Trt + Chr_Ccc.Prt * Chr_Ccc.Rrt + Chr_Ccc.Pat 
               * Chr_Ccc.Art + Chr_Ccc.Ppt * Chr_Ccc.Prt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Trt + Chr_Ccc.Prt * Chr_Ccc.Rrt + Chr_Ccc.Pat 
               * Chr_Ccc.Art + Chr_Ccc.Ppt * Chr_Ccc.Prt );
 R_Cccc.Prtr = _d_Chr_Cccc.Prrt - _d_Chr_Cccc.Prtr + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Trr + Chr_Ccc.Prt * Chr_Ccc.Rrr + Chr_Ccc.Pat 
               * Chr_Ccc.Arr + Chr_Ccc.Ppt * Chr_Ccc.Prr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Trt + Chr_Ccc.Prr * Chr_Ccc.Rrt + Chr_Ccc.Par 
               * Chr_Ccc.Art + Chr_Ccc.Ppr * Chr_Ccc.Prt );
 R_Cccc.Prta = _d_Chr_Cccc.Prat - _d_Chr_Cccc.Prta + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tra + Chr_Ccc.Prt * Chr_Ccc.Rra + Chr_Ccc.Pat 
               * Chr_Ccc.Ara + Chr_Ccc.Ppt * Chr_Ccc.Pra ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Trt + Chr_Ccc.Pra * Chr_Ccc.Rrt + Chr_Ccc.Paa 
               * Chr_Ccc.Art + Chr_Ccc.Ppa * Chr_Ccc.Prt );
 R_Cccc.Prtp = _d_Chr_Cccc.Prpt - _d_Chr_Cccc.Prtp + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Trp + Chr_Ccc.Prt * Chr_Ccc.Rrp + Chr_Ccc.Pat 
               * Chr_Ccc.Arp + Chr_Ccc.Ppt * Chr_Ccc.Prp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Trt + Chr_Ccc.Prp * Chr_Ccc.Rrt + Chr_Ccc.Pap 
               * Chr_Ccc.Art + Chr_Ccc.Ppp * Chr_Ccc.Prt );
 R_Cccc.Prrt = _d_Chr_Cccc.Prtr - _d_Chr_Cccc.Prrt + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Trt + Chr_Ccc.Prr * Chr_Ccc.Rrt + Chr_Ccc.Par 
               * Chr_Ccc.Art + Chr_Ccc.Ppr * Chr_Ccc.Prt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Trr + Chr_Ccc.Prt * Chr_Ccc.Rrr + Chr_Ccc.Pat 
               * Chr_Ccc.Arr + Chr_Ccc.Ppt * Chr_Ccc.Prr );
 R_Cccc.Prrr = _d_Chr_Cccc.Prrr - _d_Chr_Cccc.Prrr + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Trr + Chr_Ccc.Prr * Chr_Ccc.Rrr + Chr_Ccc.Par 
               * Chr_Ccc.Arr + Chr_Ccc.Ppr * Chr_Ccc.Prr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Trr + Chr_Ccc.Prr * Chr_Ccc.Rrr + Chr_Ccc.Par 
               * Chr_Ccc.Arr + Chr_Ccc.Ppr * Chr_Ccc.Prr );
 R_Cccc.Prra = _d_Chr_Cccc.Prar - _d_Chr_Cccc.Prra + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tra + Chr_Ccc.Prr * Chr_Ccc.Rra + Chr_Ccc.Par 
               * Chr_Ccc.Ara + Chr_Ccc.Ppr * Chr_Ccc.Pra ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Trr + Chr_Ccc.Pra * Chr_Ccc.Rrr + Chr_Ccc.Paa 
               * Chr_Ccc.Arr + Chr_Ccc.Ppa * Chr_Ccc.Prr );
 R_Cccc.Prrp = _d_Chr_Cccc.Prpr - _d_Chr_Cccc.Prrp + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Trp + Chr_Ccc.Prr * Chr_Ccc.Rrp + Chr_Ccc.Par 
               * Chr_Ccc.Arp + Chr_Ccc.Ppr * Chr_Ccc.Prp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Trr + Chr_Ccc.Prp * Chr_Ccc.Rrr + Chr_Ccc.Pap 
               * Chr_Ccc.Arr + Chr_Ccc.Ppp * Chr_Ccc.Prr );
 R_Cccc.Prat = _d_Chr_Cccc.Prta - _d_Chr_Cccc.Prat + ( Chr_Ccc.Pta 
               * Chr_Ccc.Trt + Chr_Ccc.Pra * Chr_Ccc.Rrt + Chr_Ccc.Paa 
               * Chr_Ccc.Art + Chr_Ccc.Ppa * Chr_Ccc.Prt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tra + Chr_Ccc.Prt * Chr_Ccc.Rra + Chr_Ccc.Pat 
               * Chr_Ccc.Ara + Chr_Ccc.Ppt * Chr_Ccc.Pra );
 R_Cccc.Prar = _d_Chr_Cccc.Prra - _d_Chr_Cccc.Prar + ( Chr_Ccc.Pta 
               * Chr_Ccc.Trr + Chr_Ccc.Pra * Chr_Ccc.Rrr + Chr_Ccc.Paa 
               * Chr_Ccc.Arr + Chr_Ccc.Ppa * Chr_Ccc.Prr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tra + Chr_Ccc.Prr * Chr_Ccc.Rra + Chr_Ccc.Par 
               * Chr_Ccc.Ara + Chr_Ccc.Ppr * Chr_Ccc.Pra );
 R_Cccc.Praa = _d_Chr_Cccc.Praa - _d_Chr_Cccc.Praa + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tra + Chr_Ccc.Pra * Chr_Ccc.Rra + Chr_Ccc.Paa 
               * Chr_Ccc.Ara + Chr_Ccc.Ppa * Chr_Ccc.Pra ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tra + Chr_Ccc.Pra * Chr_Ccc.Rra + Chr_Ccc.Paa 
               * Chr_Ccc.Ara + Chr_Ccc.Ppa * Chr_Ccc.Pra );
 R_Cccc.Prap = _d_Chr_Cccc.Prpa - _d_Chr_Cccc.Prap + ( Chr_Ccc.Pta 
               * Chr_Ccc.Trp + Chr_Ccc.Pra * Chr_Ccc.Rrp + Chr_Ccc.Paa 
               * Chr_Ccc.Arp + Chr_Ccc.Ppa * Chr_Ccc.Prp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tra + Chr_Ccc.Prp * Chr_Ccc.Rra + Chr_Ccc.Pap 
               * Chr_Ccc.Ara + Chr_Ccc.Ppp * Chr_Ccc.Pra );
 R_Cccc.Prpt = _d_Chr_Cccc.Prtp - _d_Chr_Cccc.Prpt + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Trt + Chr_Ccc.Prp * Chr_Ccc.Rrt + Chr_Ccc.Pap 
               * Chr_Ccc.Art + Chr_Ccc.Ppp * Chr_Ccc.Prt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Trp + Chr_Ccc.Prt * Chr_Ccc.Rrp + Chr_Ccc.Pat 
               * Chr_Ccc.Arp + Chr_Ccc.Ppt * Chr_Ccc.Prp );
 R_Cccc.Prpr = _d_Chr_Cccc.Prrp - _d_Chr_Cccc.Prpr + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Trr + Chr_Ccc.Prp * Chr_Ccc.Rrr + Chr_Ccc.Pap 
               * Chr_Ccc.Arr + Chr_Ccc.Ppp * Chr_Ccc.Prr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Trp + Chr_Ccc.Prr * Chr_Ccc.Rrp + Chr_Ccc.Par 
               * Chr_Ccc.Arp + Chr_Ccc.Ppr * Chr_Ccc.Prp );
 R_Cccc.Prpa = _d_Chr_Cccc.Prap - _d_Chr_Cccc.Prpa + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tra + Chr_Ccc.Prp * Chr_Ccc.Rra + Chr_Ccc.Pap 
               * Chr_Ccc.Ara + Chr_Ccc.Ppp * Chr_Ccc.Pra ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Trp + Chr_Ccc.Pra * Chr_Ccc.Rrp + Chr_Ccc.Paa 
               * Chr_Ccc.Arp + Chr_Ccc.Ppa * Chr_Ccc.Prp );
 R_Cccc.Prpp = _d_Chr_Cccc.Prpp - _d_Chr_Cccc.Prpp + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Trp + Chr_Ccc.Prp * Chr_Ccc.Rrp + Chr_Ccc.Pap 
               * Chr_Ccc.Arp + Chr_Ccc.Ppp * Chr_Ccc.Prp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Trp + Chr_Ccc.Prp * Chr_Ccc.Rrp + Chr_Ccc.Pap 
               * Chr_Ccc.Arp + Chr_Ccc.Ppp * Chr_Ccc.Prp );
 R_Cccc.Patt = _d_Chr_Cccc.Patt - _d_Chr_Cccc.Patt + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tat + Chr_Ccc.Prt * Chr_Ccc.Rat + Chr_Ccc.Pat 
               * Chr_Ccc.Aat + Chr_Ccc.Ppt * Chr_Ccc.Pat ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tat + Chr_Ccc.Prt * Chr_Ccc.Rat + Chr_Ccc.Pat 
               * Chr_Ccc.Aat + Chr_Ccc.Ppt * Chr_Ccc.Pat );
 R_Cccc.Patr = _d_Chr_Cccc.Part - _d_Chr_Cccc.Patr + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tar + Chr_Ccc.Prt * Chr_Ccc.Rar + Chr_Ccc.Pat 
               * Chr_Ccc.Aar + Chr_Ccc.Ppt * Chr_Ccc.Par ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tat + Chr_Ccc.Prr * Chr_Ccc.Rat + Chr_Ccc.Par 
               * Chr_Ccc.Aat + Chr_Ccc.Ppr * Chr_Ccc.Pat );
 R_Cccc.Pata = _d_Chr_Cccc.Paat - _d_Chr_Cccc.Pata + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Taa + Chr_Ccc.Prt * Chr_Ccc.Raa + Chr_Ccc.Pat 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppt * Chr_Ccc.Paa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tat + Chr_Ccc.Pra * Chr_Ccc.Rat + Chr_Ccc.Paa 
               * Chr_Ccc.Aat + Chr_Ccc.Ppa * Chr_Ccc.Pat );
 R_Cccc.Patp = _d_Chr_Cccc.Papt - _d_Chr_Cccc.Patp + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tap + Chr_Ccc.Prt * Chr_Ccc.Rap + Chr_Ccc.Pat 
               * Chr_Ccc.Aap + Chr_Ccc.Ppt * Chr_Ccc.Pap ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tat + Chr_Ccc.Prp * Chr_Ccc.Rat + Chr_Ccc.Pap 
               * Chr_Ccc.Aat + Chr_Ccc.Ppp * Chr_Ccc.Pat );
 R_Cccc.Part = _d_Chr_Cccc.Patr - _d_Chr_Cccc.Part + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tat + Chr_Ccc.Prr * Chr_Ccc.Rat + Chr_Ccc.Par 
               * Chr_Ccc.Aat + Chr_Ccc.Ppr * Chr_Ccc.Pat ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tar + Chr_Ccc.Prt * Chr_Ccc.Rar + Chr_Ccc.Pat 
               * Chr_Ccc.Aar + Chr_Ccc.Ppt * Chr_Ccc.Par );
 R_Cccc.Parr = _d_Chr_Cccc.Parr - _d_Chr_Cccc.Parr + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tar + Chr_Ccc.Prr * Chr_Ccc.Rar + Chr_Ccc.Par 
               * Chr_Ccc.Aar + Chr_Ccc.Ppr * Chr_Ccc.Par ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tar + Chr_Ccc.Prr * Chr_Ccc.Rar + Chr_Ccc.Par 
               * Chr_Ccc.Aar + Chr_Ccc.Ppr * Chr_Ccc.Par );
 R_Cccc.Para = _d_Chr_Cccc.Paar - _d_Chr_Cccc.Para + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Taa + Chr_Ccc.Prr * Chr_Ccc.Raa + Chr_Ccc.Par 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppr * Chr_Ccc.Paa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tar + Chr_Ccc.Pra * Chr_Ccc.Rar + Chr_Ccc.Paa 
               * Chr_Ccc.Aar + Chr_Ccc.Ppa * Chr_Ccc.Par );
 R_Cccc.Parp = _d_Chr_Cccc.Papr - _d_Chr_Cccc.Parp + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tap + Chr_Ccc.Prr * Chr_Ccc.Rap + Chr_Ccc.Par 
               * Chr_Ccc.Aap + Chr_Ccc.Ppr * Chr_Ccc.Pap ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tar + Chr_Ccc.Prp * Chr_Ccc.Rar + Chr_Ccc.Pap 
               * Chr_Ccc.Aar + Chr_Ccc.Ppp * Chr_Ccc.Par );
 R_Cccc.Paat = _d_Chr_Cccc.Pata - _d_Chr_Cccc.Paat + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tat + Chr_Ccc.Pra * Chr_Ccc.Rat + Chr_Ccc.Paa 
               * Chr_Ccc.Aat + Chr_Ccc.Ppa * Chr_Ccc.Pat ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Taa + Chr_Ccc.Prt * Chr_Ccc.Raa + Chr_Ccc.Pat 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppt * Chr_Ccc.Paa );
 R_Cccc.Paar = _d_Chr_Cccc.Para - _d_Chr_Cccc.Paar + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tar + Chr_Ccc.Pra * Chr_Ccc.Rar + Chr_Ccc.Paa 
               * Chr_Ccc.Aar + Chr_Ccc.Ppa * Chr_Ccc.Par ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Taa + Chr_Ccc.Prr * Chr_Ccc.Raa + Chr_Ccc.Par 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppr * Chr_Ccc.Paa );
 R_Cccc.Paaa = _d_Chr_Cccc.Paaa - _d_Chr_Cccc.Paaa + ( Chr_Ccc.Pta 
               * Chr_Ccc.Taa + Chr_Ccc.Pra * Chr_Ccc.Raa + Chr_Ccc.Paa 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppa * Chr_Ccc.Paa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Taa + Chr_Ccc.Pra * Chr_Ccc.Raa + Chr_Ccc.Paa 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppa * Chr_Ccc.Paa );
 R_Cccc.Paap = _d_Chr_Cccc.Papa - _d_Chr_Cccc.Paap + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tap + Chr_Ccc.Pra * Chr_Ccc.Rap + Chr_Ccc.Paa 
               * Chr_Ccc.Aap + Chr_Ccc.Ppa * Chr_Ccc.Pap ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Taa + Chr_Ccc.Prp * Chr_Ccc.Raa + Chr_Ccc.Pap 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppp * Chr_Ccc.Paa );
 R_Cccc.Papt = _d_Chr_Cccc.Patp - _d_Chr_Cccc.Papt + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tat + Chr_Ccc.Prp * Chr_Ccc.Rat + Chr_Ccc.Pap 
               * Chr_Ccc.Aat + Chr_Ccc.Ppp * Chr_Ccc.Pat ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tap + Chr_Ccc.Prt * Chr_Ccc.Rap + Chr_Ccc.Pat 
               * Chr_Ccc.Aap + Chr_Ccc.Ppt * Chr_Ccc.Pap );
 R_Cccc.Papr = _d_Chr_Cccc.Parp - _d_Chr_Cccc.Papr + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tar + Chr_Ccc.Prp * Chr_Ccc.Rar + Chr_Ccc.Pap 
               * Chr_Ccc.Aar + Chr_Ccc.Ppp * Chr_Ccc.Par ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tap + Chr_Ccc.Prr * Chr_Ccc.Rap + Chr_Ccc.Par 
               * Chr_Ccc.Aap + Chr_Ccc.Ppr * Chr_Ccc.Pap );
 R_Cccc.Papa = _d_Chr_Cccc.Paap - _d_Chr_Cccc.Papa + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Taa + Chr_Ccc.Prp * Chr_Ccc.Raa + Chr_Ccc.Pap 
               * Chr_Ccc.Aaa + Chr_Ccc.Ppp * Chr_Ccc.Paa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tap + Chr_Ccc.Pra * Chr_Ccc.Rap + Chr_Ccc.Paa 
               * Chr_Ccc.Aap + Chr_Ccc.Ppa * Chr_Ccc.Pap );
 R_Cccc.Papp = _d_Chr_Cccc.Papp - _d_Chr_Cccc.Papp + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tap + Chr_Ccc.Prp * Chr_Ccc.Rap + Chr_Ccc.Pap 
               * Chr_Ccc.Aap + Chr_Ccc.Ppp * Chr_Ccc.Pap ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tap + Chr_Ccc.Prp * Chr_Ccc.Rap + Chr_Ccc.Pap 
               * Chr_Ccc.Aap + Chr_Ccc.Ppp * Chr_Ccc.Pap );
 R_Cccc.Pptt = _d_Chr_Cccc.Pptt - _d_Chr_Cccc.Pptt + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpt + Chr_Ccc.Prt * Chr_Ccc.Rpt + Chr_Ccc.Pat 
               * Chr_Ccc.Apt + Chr_Ccc.Ppt * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpt + Chr_Ccc.Prt * Chr_Ccc.Rpt + Chr_Ccc.Pat 
               * Chr_Ccc.Apt + Chr_Ccc.Ppt * Chr_Ccc.Ppt );
 R_Cccc.Pptr = _d_Chr_Cccc.Pprt - _d_Chr_Cccc.Pptr + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpr + Chr_Ccc.Prt * Chr_Ccc.Rpr + Chr_Ccc.Pat 
               * Chr_Ccc.Apr + Chr_Ccc.Ppt * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpt + Chr_Ccc.Prr * Chr_Ccc.Rpt + Chr_Ccc.Par 
               * Chr_Ccc.Apt + Chr_Ccc.Ppr * Chr_Ccc.Ppt );
 R_Cccc.Ppta = _d_Chr_Cccc.Ppat - _d_Chr_Cccc.Ppta + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpa + Chr_Ccc.Prt * Chr_Ccc.Rpa + Chr_Ccc.Pat 
               * Chr_Ccc.Apa + Chr_Ccc.Ppt * Chr_Ccc.Ppa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpt + Chr_Ccc.Pra * Chr_Ccc.Rpt + Chr_Ccc.Paa 
               * Chr_Ccc.Apt + Chr_Ccc.Ppa * Chr_Ccc.Ppt );
 R_Cccc.Pptp = _d_Chr_Cccc.Pppt - _d_Chr_Cccc.Pptp + ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpp + Chr_Ccc.Prt * Chr_Ccc.Rpp + Chr_Ccc.Pat 
               * Chr_Ccc.App + Chr_Ccc.Ppt * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpt + Chr_Ccc.Prp * Chr_Ccc.Rpt + Chr_Ccc.Pap 
               * Chr_Ccc.Apt + Chr_Ccc.Ppp * Chr_Ccc.Ppt );
 R_Cccc.Pprt = _d_Chr_Cccc.Pptr - _d_Chr_Cccc.Pprt + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpt + Chr_Ccc.Prr * Chr_Ccc.Rpt + Chr_Ccc.Par 
               * Chr_Ccc.Apt + Chr_Ccc.Ppr * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpr + Chr_Ccc.Prt * Chr_Ccc.Rpr + Chr_Ccc.Pat 
               * Chr_Ccc.Apr + Chr_Ccc.Ppt * Chr_Ccc.Ppr );
 R_Cccc.Pprr = _d_Chr_Cccc.Pprr - _d_Chr_Cccc.Pprr + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpr + Chr_Ccc.Prr * Chr_Ccc.Rpr + Chr_Ccc.Par 
               * Chr_Ccc.Apr + Chr_Ccc.Ppr * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpr + Chr_Ccc.Prr * Chr_Ccc.Rpr + Chr_Ccc.Par 
               * Chr_Ccc.Apr + Chr_Ccc.Ppr * Chr_Ccc.Ppr );
 R_Cccc.Ppra = _d_Chr_Cccc.Ppar - _d_Chr_Cccc.Ppra + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpa + Chr_Ccc.Prr * Chr_Ccc.Rpa + Chr_Ccc.Par 
               * Chr_Ccc.Apa + Chr_Ccc.Ppr * Chr_Ccc.Ppa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpr + Chr_Ccc.Pra * Chr_Ccc.Rpr + Chr_Ccc.Paa 
               * Chr_Ccc.Apr + Chr_Ccc.Ppa * Chr_Ccc.Ppr );
 R_Cccc.Pprp = _d_Chr_Cccc.Pppr - _d_Chr_Cccc.Pprp + ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpp + Chr_Ccc.Prr * Chr_Ccc.Rpp + Chr_Ccc.Par 
               * Chr_Ccc.App + Chr_Ccc.Ppr * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpr + Chr_Ccc.Prp * Chr_Ccc.Rpr + Chr_Ccc.Pap 
               * Chr_Ccc.Apr + Chr_Ccc.Ppp * Chr_Ccc.Ppr );
 R_Cccc.Ppat = _d_Chr_Cccc.Ppta - _d_Chr_Cccc.Ppat + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpt + Chr_Ccc.Pra * Chr_Ccc.Rpt + Chr_Ccc.Paa 
               * Chr_Ccc.Apt + Chr_Ccc.Ppa * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpa + Chr_Ccc.Prt * Chr_Ccc.Rpa + Chr_Ccc.Pat 
               * Chr_Ccc.Apa + Chr_Ccc.Ppt * Chr_Ccc.Ppa );
 R_Cccc.Ppar = _d_Chr_Cccc.Ppra - _d_Chr_Cccc.Ppar + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpr + Chr_Ccc.Pra * Chr_Ccc.Rpr + Chr_Ccc.Paa 
               * Chr_Ccc.Apr + Chr_Ccc.Ppa * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpa + Chr_Ccc.Prr * Chr_Ccc.Rpa + Chr_Ccc.Par 
               * Chr_Ccc.Apa + Chr_Ccc.Ppr * Chr_Ccc.Ppa );
 R_Cccc.Ppaa = _d_Chr_Cccc.Ppaa - _d_Chr_Cccc.Ppaa + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpa + Chr_Ccc.Pra * Chr_Ccc.Rpa + Chr_Ccc.Paa 
               * Chr_Ccc.Apa + Chr_Ccc.Ppa * Chr_Ccc.Ppa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpa + Chr_Ccc.Pra * Chr_Ccc.Rpa + Chr_Ccc.Paa 
               * Chr_Ccc.Apa + Chr_Ccc.Ppa * Chr_Ccc.Ppa );
 R_Cccc.Ppap = _d_Chr_Cccc.Pppa - _d_Chr_Cccc.Ppap + ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpp + Chr_Ccc.Pra * Chr_Ccc.Rpp + Chr_Ccc.Paa 
               * Chr_Ccc.App + Chr_Ccc.Ppa * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpa + Chr_Ccc.Prp * Chr_Ccc.Rpa + Chr_Ccc.Pap 
               * Chr_Ccc.Apa + Chr_Ccc.Ppp * Chr_Ccc.Ppa );
 R_Cccc.Pppt = _d_Chr_Cccc.Pptp - _d_Chr_Cccc.Pppt + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpt + Chr_Ccc.Prp * Chr_Ccc.Rpt + Chr_Ccc.Pap 
               * Chr_Ccc.Apt + Chr_Ccc.Ppp * Chr_Ccc.Ppt ) - ( Chr_Ccc.Ptt 
               * Chr_Ccc.Tpp + Chr_Ccc.Prt * Chr_Ccc.Rpp + Chr_Ccc.Pat 
               * Chr_Ccc.App + Chr_Ccc.Ppt * Chr_Ccc.Ppp );
 R_Cccc.Pppr = _d_Chr_Cccc.Pprp - _d_Chr_Cccc.Pppr + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpr + Chr_Ccc.Prp * Chr_Ccc.Rpr + Chr_Ccc.Pap 
               * Chr_Ccc.Apr + Chr_Ccc.Ppp * Chr_Ccc.Ppr ) - ( Chr_Ccc.Ptr 
               * Chr_Ccc.Tpp + Chr_Ccc.Prr * Chr_Ccc.Rpp + Chr_Ccc.Par 
               * Chr_Ccc.App + Chr_Ccc.Ppr * Chr_Ccc.Ppp );
 R_Cccc.Pppa = _d_Chr_Cccc.Ppap - _d_Chr_Cccc.Pppa + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpa + Chr_Ccc.Prp * Chr_Ccc.Rpa + Chr_Ccc.Pap 
               * Chr_Ccc.Apa + Chr_Ccc.Ppp * Chr_Ccc.Ppa ) - ( Chr_Ccc.Pta 
               * Chr_Ccc.Tpp + Chr_Ccc.Pra * Chr_Ccc.Rpp + Chr_Ccc.Paa 
               * Chr_Ccc.App + Chr_Ccc.Ppa * Chr_Ccc.Ppp );
 R_Cccc.Pppp = _d_Chr_Cccc.Pppp - _d_Chr_Cccc.Pppp + ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpp + Chr_Ccc.Prp * Chr_Ccc.Rpp + Chr_Ccc.Pap 
               * Chr_Ccc.App + Chr_Ccc.Ppp * Chr_Ccc.Ppp ) - ( Chr_Ccc.Ptp 
               * Chr_Ccc.Tpp + Chr_Ccc.Prp * Chr_Ccc.Rpp + Chr_Ccc.Pap 
               * Chr_Ccc.App + Chr_Ccc.Ppp * Chr_Ccc.Ppp );

 R_cc.tt = ( R_Cccc.Tttt + R_Cccc.Rtrt + R_Cccc.Atat + R_Cccc.Ptpt );
 R_cc.tr = ( R_Cccc.Tttr + R_Cccc.Rtrr + R_Cccc.Atar + R_Cccc.Ptpr );
 R_cc.ta = ( R_Cccc.Ttta + R_Cccc.Rtra + R_Cccc.Ataa + R_Cccc.Ptpa );
 R_cc.tp = ( R_Cccc.Tttp + R_Cccc.Rtrp + R_Cccc.Atap + R_Cccc.Ptpp );
 R_cc.rt = ( R_Cccc.Trtt + R_Cccc.Rrrt + R_Cccc.Arat + R_Cccc.Prpt );
 R_cc.rr = ( R_Cccc.Trtr + R_Cccc.Rrrr + R_Cccc.Arar + R_Cccc.Prpr );
 R_cc.ra = ( R_Cccc.Trta + R_Cccc.Rrra + R_Cccc.Araa + R_Cccc.Prpa );
 R_cc.rp = ( R_Cccc.Trtp + R_Cccc.Rrrp + R_Cccc.Arap + R_Cccc.Prpp );
 R_cc.at = ( R_Cccc.Tatt + R_Cccc.Rart + R_Cccc.Aaat + R_Cccc.Papt );
 R_cc.ar = ( R_Cccc.Tatr + R_Cccc.Rarr + R_Cccc.Aaar + R_Cccc.Papr );
 R_cc.aa = ( R_Cccc.Tata + R_Cccc.Rara + R_Cccc.Aaaa + R_Cccc.Papa );
 R_cc.ap = ( R_Cccc.Tatp + R_Cccc.Rarp + R_Cccc.Aaap + R_Cccc.Papp );
 R_cc.pt = ( R_Cccc.Tptt + R_Cccc.Rprt + R_Cccc.Apat + R_Cccc.Pppt );
 R_cc.pr = ( R_Cccc.Tptr + R_Cccc.Rprr + R_Cccc.Apar + R_Cccc.Pppr );
 R_cc.pa = ( R_Cccc.Tpta + R_Cccc.Rpra + R_Cccc.Apaa + R_Cccc.Pppa );
 R_cc.pp = ( R_Cccc.Tptp + R_Cccc.Rprp + R_Cccc.Apap + R_Cccc.Pppp );

 R = ( g_CC.TT * R_cc.tt + g_CC.TR * R_cc.rt + g_CC.TA * R_cc.at + g_CC.TP 
     * R_cc.pt + g_CC.RT * R_cc.tr + g_CC.RR * R_cc.rr + g_CC.RA * R_cc.ar 
     + g_CC.RP * R_cc.pr + g_CC.AT * R_cc.ta + g_CC.AR * R_cc.ra + g_CC.AA 
     * R_cc.aa + g_CC.AP * R_cc.pa + g_CC.PT * R_cc.tp + g_CC.PR * R_cc.rp 
     + g_CC.PA * R_cc.ap + g_CC.PP * R_cc.pp );

/* Close GRPP Block */
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


