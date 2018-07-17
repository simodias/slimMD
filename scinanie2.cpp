//ostatnia modyfikacja 02. 03. 2013

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <conio.h>
#include <cstdlib>
#include <stdlib.h>
#include <string>
#include <direct.h>


//------------------------------------------------------------------------------
//FUNKCJE MATEMATYCZNE
//------------------------------------------------------------------------------

//dzielenie---------------------------------------------------------------------
double dziel (double D1, double D2)
{
double iloraz;
iloraz=D1*pow(D2,-1);
return iloraz;       
}
//pierwiastkowanie--------------------------------------------------------------
double pierw (double D1)
{
double pierwiastek;
pierwiastek=pow(D1,0.5);
return pierwiastek;       
}
//kwadratowanie-----------------------------------------------------------------
double kw (double D1)
{
double kwadrat;
kwadrat=pow(D1,2);
return kwadrat;       
}
//zaokraglanie------------------------------------------------------------------
double przybliz (double D1)
{
double zaokraglenie;
zaokraglenie=int(D1+0.5);
return zaokraglenie;
} 

//------------------------------------------------------------------------------
//DEKLARACJE ZMIENNYCH
//------------------------------------------------------------------------------
#include "stale.h"

int NATIN;
int READ, wczytaj_sciany, granica_odczytu_scian;
double RX[NATOMS], RY[NATOMS], RZ[NATOMS], RXR[NATOMS], RYR[NATOMS], RZR[NATOMS];
double DELTARX, DELTARY, DELTARZ, RX_old[NATOMS], RY_old[NATOMS], RZ_old[NATOMS];
double RX_p[NATOMS], RY_p[NATOMS], RZ_p[NATOMS], RX_eq[NATOMS], RY_eq[NATOMS], RZ_eq[NATOMS];

double VX[NATOMS], VY[NATOMS], VZ[NATOMS], VX_old[NATOMS], VY_old[NATOMS], VZ_old[NATOMS];
double VX_p[NATOMS], VY_p[NATOMS], VZ_p[NATOMS], VX_pec[NATOMS], VY_pec[NATOMS], VZ_pec[NATOMS], VPEC[NATOMS];

double FX[NATOMS], FY[NATOMS], FZ[NATOMS], FX_old[NATOMS], FY_old[NATOMS], FZ_old[NATOMS];
double FX_p[NATOMS], FY_p[NATOMS], FZ_p[NATOMS], DFX[NATOMS], DFY[NATOMS], DFZ[NATOMS], FTET[NATOMS];

double U[NATOMS], DU[NATOMS];

double VXH2[NATOMS], VXH1[NATOMS], VYH2[NATOMS], VYH1[NATOMS], VZH2[NATOMS], VZH1[NATOMS];

double KB, T0, T0_down, T0_up, DT, T, T_down, T_up, T_mid, T_mid_2, T_walls, Q, Q_in, tim;
double X, Y, Z, RCUT, WCA, R, R2, RI, a, b, EUP, FXP, FYP, FZP, EPC, EPC_walls, EPC_inter, EPC_w_down;
double EPC_w_up, sumEPC_inter, sumEPC_w_down, sumEPC_w_up, meanEPC_inter, meanEPC_w_down;
double meanEPC_w_up, EKC_w_down, EKC_w_up, EKC, EKC_walls, EKC_inter, K, shift, N1, N2; 
double calka, calka_walls, calka_w_down, calka_w_up, calka_inter, TETF, TETF2;
double DST, SP1, SP2, DBOX, DBOZ, GRAN, rho, Volume, TOTPX, TOTPY, TOTPZ;
double EnergyQ, EnergyQ2, EnergyU, EnergyK, Energy;
double DX, DY, DZ, DR;
double Pt, P0, Qp, vzz, vzz_dec, F_down, F_up, V_down, V_up, sumP, k4, k6;
double AX, AY, AZ, AR;
double PXX, PYY, PZZ, PXY, PYX, PXZ, PZX, PZY, PYZ, PVIR;
double sumPXX, sumPYY, sumPZZ, sumPXY, sumPYX, sumPXZ, sumPZX, sumPZY, sumPYZ;
double meanPXX, meanPYY, meanPZZ, meanPXY, meanPYX, meanPXZ, meanPZX, meanPZY, meanPYZ;
double FUP2, FDOWN2;
int TERM, XPERIOD, YPERIOD, ZPERIOD, BAR;
int NUNIT, M, g, g_in, RUN, LOOP, JUMP, JUMP2;
double RUN_TIME;
int wskA, wskB, wskC, wskD, wskE, wskF, wskG, wskH, SURFSTART, SURFEND, LEYEND, ITSURF, INTI, INTJ;
int counter, marker[NATOMS], NLL, list[NATOMS][NATOMS], countX, ccc;
int StartVelocities;
double SKIN;
double DRX[NATOMS], DRY[NATOMS], DRZ[NATOMS]; 
double SUM_EP, SUM_T, SUM_T_up, SUM_T_down, SUM_T_mid, SUM_T_mid_2, SUM_T_walls, SUM_P; 
double meanEP, meanT, meanT_up, meanT_down, meanT_mid, meanT_mid_2, meanT_walls, meanP, meanPZ, DENOM;
double scale;
double CPX, CPY, CPZ, CMPX, CMPY, CMPZ;
double tau;
int cleaning;
double FXRIGHT, FXLEFT, FYRIGHT, FYLEFT, FZUP, FZDOWN, PTX, PTY, PTZ, PVIRX, PVIRY, PVIRZ;
double sumPTX, sumPTY, sumPTZ;
int LX[NATOMS], LY[NATOMS];
double DSTX, DSTY, DSTZ, BX, BY, BZ, BZ_old, GAP, CENTGAP, DIMZ;
int LEY, CONSTR;
double PVZ, EKZ;
double sumPVZ, sumDIM;
double VSH;
int SHEAR;
int DIFINT, DIFCON, NSURF;
double DIMRATE;
double SURFX_UP[2*NUNITX*NUNITY], SURFY_UP[2*NUNITX*NUNITY], SURFZ_UP[2*NUNITX*NUNITY];
double SURFX_DOWN[2*NUNITX*NUNITY], SURFY_DOWN[2*NUNITX*NUNITY], SURFZ_DOWN[2*NUNITX*NUNITY];
double mD, mU, RaD, RaU, sigma2D, sigma2U, SkD, SkU, KurtD, KurtU;
double fric1, fric2, FBXU, FBXD, FBYU, FBYD;
double vzz_old;
double E, F, A, MA, MB, C;
double tauH;
double phi[NATOMS], CR, CALFA, SUMCR, SUMCRCA;
double dx, dy, dz, rpp, r0, r1, ys;
double counterCR1, counterCR2, counterCR3, counterCRCA;
double MI1, sumMI1, meanMI1, MI2, sumMI2, meanMI2, MI3, sumMI3, meanMI3, MI4, sumMI4, meanMI4, MI5, sumMI5, meanMI5, MI6, sumMI6, meanMI6, SUMFRIC1, SUMFRIC2, SUMFRIC3, SUMFRIC4, SUMFRIC5, SUMFRIC6;
double POSZ1, POSZ2, LENGTH, sumL, meanL;
int iterj, bin, WHZ, GET, binHV;
int histX[granica], histY[granica], histZ[granica], histZ2[granica], NhistZ[granica], BAR_CHOICE;
int termX, termY, termZ, termWall, odczyt_wejscie;
double mass;
double calkaBEGIN, DH, meanDH, sumDH;
double CPXA, CPYA, CPZA, CPXB, CPYB, CPZB, CPXC, CPYC, CPZC;
double SMX, SMY, SMZ, SMXA, SMYA, SMZA, SMXB, SMYB, SMZB, SMXC, SMYC, SMZC;
double sumSMX, sumSMY, sumSMZ, sumSMXA, sumSMYA, sumSMZA, sumSMXB, sumSMYB, sumSMZB, sumSMXC, sumSMYC, sumSMZC;
double meanSMX, meanSMY, meanSMZ, meanSMXA, meanSMYA, meanSMZA, meanSMXB, meanSMYB, meanSMZB, meanSMXC, meanSMYC, meanSMZC;
double UAA, UA, UBB, UAB, UBC, UCC, UC, ULJ[NATOMS], UPH[NATOMS], ELJ, EPH;
double FAAX, FAX, FBBX, FABX, FBCX, FCCX, FCX;
double FAAY, FAY, FBBY, FABY, FBCY, FCCY, FCY;
double FAAZ, FAZ, FBBZ, FABZ1, FABZ2, FBCZ1, FBCZ2, FCCZ, FCZ;
double cax, cay, caz;
double crc;
double faktor_k;
double SUM_DR, meanDR;
double SXP, SYP, SZP;
double PDIM1, PDIM2, DIMZ2, DIMZ3, Volum1, Volum2, Volum3, Volum4, sumVolum1, sumVolum2, sumVolum3, meanVolum1, meanVolum2, meanVolum3;
double VIR, PVIR1, PVIR2, PVIR3, PVIR4, SPVIR1, SPVIR2, SPVIR3, NEWP1, NEWP2; 
double sumPVIR1, sumPVIR2, sumPVIR3, sumPVIR4, sumSPVIR1, sumSPVIR2, sumSPVIR3, sumNEWP1, sumNEWP2; 
double meanPVIR1, meanPVIR2, meanPVIR3, meanPVIR4, meanSPVIR1, meanSPVIR2, meanSPVIR3, meanNEWP1, meanNEWP2; 

double SUM_PXX, SUM_PYY, SUM_PZZ, SUM_PXY, SUM_PYX, SUM_PZX, SUM_PXZ, SUM_PYZ, SUM_PZY;
double SUM_FDOWN2, SUM_FUP2, SUM_FBXD, SUM_FBXU, SUM_FBYD, SUM_FBYU;
double deltaP, deltaJumpP, deltaV;
int licznikomp;
double resol;
double dEW1, dEW2, sumDEW1, sumDEW2, meanDEW1, meanDEW2;
double distA, distC;
double shFA, shFC, sumShFA, sumShFC, meanShFA, meanShFC;
int cisnienieNaturalne;       

int GR[RANG], WGR[RANG], DISTANCE; 
double NGR[RANG], NWGR[RANG], TT[RANG];   
int probny_counter, ucount1, ucount2;

int binT, binD, WH, CNT, binP, binL, binRZ;
int histT[RANG], histT_up[RANG], histT_down[RANG], histT_mid[RANG], histT_mid_2[RANG], histP[RANG2], histL[RANG2];
double L0, Pb;

double tmp;
int number[NATOMS], counter_freez;

double scale_all, scale_up, scale_down;
int freq_scale, potencjal;
int scale_vel_w;
double ENKX, ENKY, ENKZ, TEMPX, TEMPY, TEMPZ, sumTX, sumTY, sumTZ, meanTX, meanTY, meanTZ;

double RTZ[TGAP], SUM_RTZ[TGAP], MEAN_RTZ[TGAP], MEAN_RTZ2[TGAP];
double RTZ_2[TGAP], SUM_RTZ_2[TGAP], MEAN_RTZ_2[TGAP];

double RTZ_X[TGAP], SUM_RTZ_X[TGAP], MEAN_RTZ_X[TGAP], MEAN_RTZ_X2[TGAP];
double RTZ_Y[TGAP], SUM_RTZ_Y[TGAP], MEAN_RTZ_Y[TGAP], MEAN_RTZ_Y2[TGAP];
double RTZ_Z[TGAP], SUM_RTZ_Z[TGAP], MEAN_RTZ_Z[TGAP], MEAN_RTZ_Z2[TGAP];

double RVX[TGAP], SUM_RTVX[TGAP], MEAN_RTVX[TGAP], MEAN_RTVX2[TGAP];
double RVY[TGAP], SUM_RTVY[TGAP], MEAN_RTVY[TGAP], MEAN_RTVY2[TGAP];
double RVY_pec[TGAP], SUM_RTVY_pec[TGAP], MEAN_RTVY_pec[TGAP], MEAN_RTVY_pec2[TGAP];
double RVZ[TGAP], SUM_RTVZ[TGAP], MEAN_RTVZ[TGAP], MEAN_RTVZ2[TGAP];

int C_VX[TGAP], C_VY[TGAP], C_VY_pec[TGAP], C_VZ[TGAP];
int COUNT_Z[TGAP], COUNT_Y[TGAP], COUNT_X[TGAP], COUNT_2[TGAP], COUNT[TGAP];

int tabDenomRTZ[TGAP], preTabDenomRTZ[TGAP];
int tabDenomRTZ_2[TGAP], preTabDenomRTZ_2[TGAP];
int tabDenomRTZ_X[TGAP], preTabDenomRTZ_X[TGAP];
int tabDenomRTZ_Y[TGAP], preTabDenomRTZ_Y[TGAP];
int tabDenomRTZ_Z[TGAP], preTabDenomRTZ_Z[TGAP];

int tabDenomRVX[TGAP], preTabDenomRVX[TGAP];
int tabDenomRVY[TGAP], preTabDenomRVY[TGAP];
int tabDenomRVY_pec[TGAP], preTabDenomRVY_pec[TGAP];
int tabDenomRVZ[TGAP], preTabDenomRVZ[TGAP];

double histVXA[TGAP], histVYA[TGAP], histVZA[TGAP], histVA[TGAP];
double histVXB[TGAP], histVYB[TGAP], histVZB[TGAP], histVB[TGAP];
double histVXC[TGAP], histVYC[TGAP], histVZC[TGAP], histVC[TGAP];

double VAP, VAPcount, sumVAP, meanVAP;
double VAP2, VAP2count, sumVAP2, meanVAP2;
double VAP3, VAP3count, sumVAP3, meanVAP3;
double VAP4, VAP4count, sumVAP4, meanVAP4;
double VAP5, VAP5count, sumVAP5, meanVAP5;
double VAP10, VAP10count, sumVAP10, meanVAP10;
double VAP20, VAP20count, sumVAP20, meanVAP20;
int VAPstart, VAPend; 
int flaga_AP, flaga_AS;

double ROZ;
double PW, sumPW, meanPW, FYext, FY0, PWX, PWY, sumPWX, sumPWY, meanPWX, meanPWY;
double VXA, VYA, VXC, VYC, sumVXA, sumVYA, sumVXC, sumVYC, meanVXA, meanVYA, meanVXC, meanVYC;
int BARO4, BARO5, SH1, SH2, SH3;
double PPXA, PPYA, PPZA;
double PPXB, PPYB, PPZB;
double PPXC, PPYC, PPZC;
double FEXT;
int ZerowaniePedu, ZerowanieSrodkaM;
double VZB5_A, VZB5_C;
double R3, KR3;
int korekta_pedu_scian;
int WR3, PAIR_TET, lista2[NATOMS][6], wsk_lista2, zerowanieSM2;
int RandA, RandB;
int dataLines;
int nis_prec, wys_prec;
double a1, a2;

double PedXY, PedYX, PedXZ, PedZX, PedZY, PedYZ;

bool swapped;
double meanPt;
double FP1U, FP1D, FP2U, FP2D, FP3U, FP3D, FP4U, FP4D, FP5U, FP5D, FP6U, FP6D; 
double FP7U, FP7D, FP8U, FP8D, FP9U, FP9D, FP10U, FP10D, FP11U, FP11D, FP12U, FP12D; 
double FP13U, FP13D, FP14U, FP14D, FP15U, FP15D, FP16U, FP16D, FP17U, FP17D, FP18U, FP18D; 
double FP19U, FP19D, FP20U, FP20D;

double FPU[PGAP], FPD[PGAP], PMOP[PGAP], sumPMOP[PGAP], meanPMOP[PGAP];

double FPU2[PGAP2], FPD2[PGAP2], PMOP2[PGAP2], PMOP_old[PGAP2], sumPMOP2[PGAP2], meanPMOP2[PGAP2], PTET[PGAP2];
double ZP;//okresla zasieg obliczania VAP
int calcGR;
double PXX_mid, PYY_mid, PZZ_mid, PYZ_mid, PVIR_mid, PVIR_mid2, PVIR_mid3;
double sumPXX_mid, sumPYY_mid, sumPZZ_mid, sumPYZ_mid;
double meanPXX_mid, meanPYY_mid, meanPZZ_mid, meanPYZ_mid;

double sumPVIR_mid1, sumPVIR_mid2, sumPVIR_mid3, meanPVIR_mid1, meanPVIR_mid2, meanPVIR_mid3; 
double T_pec, SUM_T_pec, meanT_pec;


int sygnal, mop_calc, PN_calc;
double omega;
int automat_pressing, automat_sliding;

double POZSURF[PGAP2], deltay, part_sig1, part_sig2, PKYYMOP[PGAP2], sumPKYYMOP[PGAP2], sumPKYYMOP2[PGAP2], meanPKYYMOP[PGAP2];
double daik, daic;
double rsign_i, rsign_j;
int fsign, vsign, kincount, rsign;
int count1, count2, count3, count4;
int count5, count6, count7, count8;
double FTET1, FTET2, sumFTET1, sumFTET2, meanFTET1, meanFTET2;
double FWU, FWD, VWU, VWD;
double PT1, PT2, sumPT1, sumPT2, meanPT1, meanPT2;
int n_old, n_0, sliding, pressing, jump_P0, zmiana_warunkow_T;
int sumk, liczsumk;
double kinet[PGAP2], timkinet[PGAP2],  sumPMOP3[PGAP2], sumPMOP4[PGAP2], VSH_K, P0_K, PMOP5[PGAP2], sumPMOP5[PGAP2];	
double znak;	
double t1, t2, t3, t_test; //zmienne timera i przeskalowanie ms na s, m, h
double XA, XC, VA, VC, beta;
double NEWPX1, NEWPX2, NEWPY1, NEWPY2, sumMI, MI, meanMI, ST, gam1, gam2;
double VW, sumVW, meanVW;
//Barostat RK4------------------------------------------------------------------
double P_p, P_old, g_old, g_p, KL1, KL2, KL3, KL4, KG1, KG2, KG3, KG4;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------      
//FUNKCJA GLOWNA
//------------------------------------------------------------------------------      
using namespace std;	

int main()
{
mkdir("wyniki");	
#include "nazwy_plikow.h"

//PARAMETRY SYMULACJI-----------------------------------------------------------

#include "settings.h"

nis_prec=4;
wys_prec=16;

//------------------------------------------------------------------------------
shift=-4*(pow(RCUT,-12)-pow(RCUT,-6));
if(potencjal==1)
{
shift=-4*(pow(WCA,-12)-pow(WCA,-6));
}
//------------------------------------------------------------------------------
wskA=NUNITX*NUNITY*4*LEY;
wskB=NATOMS-wskA;

wskC=wskA+1*wskA;
wskD=wskB-1*wskA;

wskE=wskA+2*wskA;//okresla przypiete atomy poddawane scinaniu
wskF=wskB-2*wskA;//okresla przypiete atomy poddawane scinaniu

wskG=wskA+3*wskA;
wskH=wskB-3*wskA;

SURFSTART=NATOMS*0.5;
SURFEND=NUNITX*NUNITY*4*0.5;
LEYEND=NUNITX*NUNITY*4;
NSURF=0.5*LEYEND;
cout<<"DR: "<<DIMRATE<<endl;
BX=pow(NATOMS*pow(DIMRATE*rho,-1),pow(3,-1));
BY=BX;
BZ=DIMRATE*BX;
BZ_old=BZ;

Volume=BX*BY*BZ;
DBOX=BX;
DBOZ=DIMRATE*BX;

DX=0; DY=0; DZ=0; DR=0; 

GET=0;
CNT=1;

//zerowanie---------------------------------------------------------------------	
flaga_AP=0;
flaga_AS=0;
wczytaj_sciany=0;
granica_odczytu_scian=0;


kincount=0;
counter_freez=1;
crc=0.0000000001;
cax=0;
cay=0;
caz=0;
sumDH=0;
probny_counter=0;
ucount1=0;
ucount2=0;
GRAN=0;
sumP=0;
SUM_EP=0;
SUM_T=0;
SUM_T_pec=0;
SUM_P=0;
DENOM=1;
tim=0;
liczsumk=0;
sumNEWP1=0;
sumNEWP2=0;
sumPT1=0;
sumPT2=0; 
a1=2.5;
a2=1;

sumVW=0;

sumTX=0;
sumTY=0;
sumTZ=0;

sumFTET1=0; 
sumFTET2=0;

sumPXX=0;
sumPYY=0;
sumPZZ=0;

sumPYZ=0;
sumPZY=0;
sumPXZ=0;
sumPZX=0;
sumPXY=0;
sumPYX=0;

sumPTX=0;
sumPTY=0;
sumPTZ=0;

sumPXX_mid=0;
sumPYY_mid=0;
sumPZZ_mid=0;
sumPYZ_mid=0;

sumDEW1=0;
sumDEW2=0;

VWU=0;
VWD=0;
VZB5_A=0;
VZB5_C=0;

sumPWX=0;
sumPWY=0;

PWX=0;
PWY=0;
distA=0; 
distC=0;

meanPWX=0;
meanPWY=0;

shFA=0;
shFC=0;
sumShFA=0;
sumShFC=0;
meanShFA=0;
meanShFC=0;

sumPVIR_mid1=0;
sumPVIR_mid2=0;
sumPVIR_mid3=0;

meanPVIR_mid1=0;
meanPVIR_mid2=0;
meanPVIR_mid3=0;

SUMFRIC1=0;
SUMFRIC2=0; 
SUMFRIC3=0;
SUMFRIC4=0;

sumPVZ=0;
sumDIM=0;
vzz_old=0;
vzz=0;
SUM_DR=0;
meanDR=0;

SUM_PXX=0;
SUM_PYY=0;
SUM_PZZ=0;
SUM_PXY=0; 
SUM_PYX=0;
SUM_PZX=0; 
SUM_PXZ=0; 
SUM_PYZ=0;
SUM_PZY=0;
SUM_FDOWN2=0;
SUM_FUP2=0;
SUM_FBXD=0;
SUM_FBXU=0;
SUM_FBYD=0;
SUM_FBYU=0;

sumEPC_inter=0;
sumEPC_w_down=0;
sumEPC_w_up=0;
meanEPC_inter=0; 
meanEPC_w_down=0;
meanEPC_w_up=0;

sumVolum1=0;
sumVolum2=0;
sumVolum3=0;

sumPVIR1=0;
sumPVIR2=0;
sumPVIR3=0;
sumPVIR4=0;
sumSPVIR1=0;
sumSPVIR2=0;
sumSPVIR3=0; 
sumVAP=0;
meanVAP=0;

sumVAP2=0;
meanVAP2=0;
sumVAP3=0;
meanVAP3=0;
sumVAP4=0;
meanVAP4=0;
sumVAP5=0;
meanVAP5=0;
sumVAP10=0;
meanVAP10=0;
sumVAP20=0;
meanVAP20=0;

sumSMX=0;
sumSMY=0;
sumSMZ=0;
sumSMXA=0;
sumSMYA=0; 
sumSMZA=0; 
sumSMXB=0;
sumSMYB=0; 
sumSMZB=0; 
sumSMXC=0; 
sumSMYC=0; 
sumSMZC=0;
sumPW=0;

sumL=0;

sumMI=0; meanMI=0; sumMI1=0; meanMI1=0; sumMI2=0; meanMI2=0;
sumMI3=0; meanMI3=0; sumMI4=0; meanMI4=0; sumMI5=0; meanMI5=0; sumMI6=0; meanMI6=0;

for(int i=0;i<NATOMS; i++)
{
	number[i]=i;
}

for (int i=0;i<TGAP;i++)
{
	RTZ[i]=0;
	COUNT[i]=0;
	SUM_RTZ[i]=0;
	MEAN_RTZ[i]=0;
	
	RTZ_2[i]=0;
	COUNT_2[i]=0;
	SUM_RTZ_2[i]=0;
	MEAN_RTZ_2[i]=0;
	
	RTZ_X[i]=0;
	COUNT_X[i]=0;
	SUM_RTZ_X[i]=0;
	MEAN_RTZ_X[i]=0;
	
	RTZ_Y[i]=0;
	COUNT_Y[i]=0;
	SUM_RTZ_Y[i]=0;
	MEAN_RTZ_Y[i]=0;
	
	RTZ_Z[i]=0;
	COUNT_Z[i]=0;
	SUM_RTZ_Z[i]=0;
	MEAN_RTZ_Z[i]=0;
	
	RVX[i]=0; 
	C_VX[i]=0;  
	SUM_RTVX[i]=0; 
	MEAN_RTVX[i]=0; 

	RVY[i]=0;  
	C_VY[i]=0;  
	SUM_RTVY[i]=0;  
	MEAN_RTVY[i]=0; 
	
	RVY_pec[i]=0;  
	C_VY_pec[i]=0;  
	SUM_RTVY_pec[i]=0;  
	MEAN_RTVY_pec[i]=0; 

	RVZ[i]=0; 
	C_VZ[i]=0; 
	SUM_RTVZ[i]=0;  
	MEAN_RTVZ[i]=0; 
	
	tabDenomRTZ[i]=0;
	preTabDenomRTZ[i]=0;
	tabDenomRTZ_2[i]=0; 
	preTabDenomRTZ_2[i]=0;
	tabDenomRTZ_X[i]=0; 
	preTabDenomRTZ_X[i]=0;
	tabDenomRTZ_Y[i]=0; 
	preTabDenomRTZ_Y[i]=0;
	tabDenomRTZ_Z[i]=0;
	preTabDenomRTZ_Z[i]=0;

	tabDenomRVX[i]=0; 
	preTabDenomRVX[i]=0;
	tabDenomRVY[i]=0;
	preTabDenomRVY[i]=0;
	tabDenomRVY_pec[i]=0; 
	preTabDenomRVY_pec[i]=0;
	tabDenomRVZ[i]=0; 
	preTabDenomRVZ[i]=0;
	
	MEAN_RTZ2[i]=0;
	MEAN_RTZ_X2[i]=0;
	MEAN_RTZ_Y2[i]=0;
	MEAN_RTZ_Z2[i]=0;
	MEAN_RTVX2[i]=0;
	MEAN_RTVY2[i]=0;
	MEAN_RTVY_pec2[i]=0;
	MEAN_RTVZ2[i]=0;
	
	histVXA[i]=0;
	histVYA[i]=0; 
	histVZA[i]=0; 
	histVA[i]=0;
	histVXB[i]=0; 
	histVYB[i]=0; 
	histVZB[i]=0; 
	histVB[i]=0;
	histVXC[i]=0; 
	histVYC[i]=0; 
	histVZC[i]=0; 
	histVC[i]=0;
	
}

for(int i=0;i<RANG;i++)
{
	GR[i]=0;
	WGR[i]=0;
	NGR[i]=0;
	NWGR[i]=0;
	histT[i]=0;
	histT_up[i]=0;
	histT_down[i]=0;
	histT_mid[i]=0;
	histT_mid_2[i]=0;
	histX[i]=0;
	histY[i]=0;
	histZ[i]=0;
	histZ2[i]=0;
}

for(int i=0;i<RANG2;i++)
{
	histP[i]=0;
	histL[i]=0;
}


for(int kk=0;kk<PGAP2;kk++)	
{
sumPMOP2[kk]=0;	
sumPMOP3[kk]=0;	
sumPMOP4[kk]=0;	
sumPMOP5[kk]=0;	
sumPKYYMOP[kk]=0;
sumPKYYMOP2[kk]=0;
PTET[kk]=0;

kinet[kk]=0;
timkinet[kk]=0;				
}
//ODCZYT PARAMETROW WEJSCIOWYCH-------------------------------------------------

double *inputData=new double[dataLines];         


if(odczyt_wejscie==1)
{

double bufor1;
fstream plikInputVk("inputVk.txt", ios::in | ios::out );      
plikInputVk >> bufor1;
VSH_K=bufor1/158;
plikInputVk.close();	
	
fstream plikInputPk("inputPk.txt", ios::in | ios::out );      
plikInputPk >> bufor1;
P0_K=bufor1/42.133;
plikInputPk.close();	
	
	
fstream plikOgolnyInput("inputGeneral.txt", ios::in | ios::out );	

	for(int i=0;i<dataLines;i++)
	{
	plikOgolnyInput >>inputData[i];
	}
	plikOgolnyInput.close();

	VSH=inputData[0]/158;
	P0=inputData[1]/42.133;
	
	DT=inputData[2];
	RUN=inputData[3];
	cleaning=inputData[4];
	a1=inputData[5];
	a2=inputData[6];
	
	BARO4=inputData[7];
	BARO5=inputData[8];
	beta=inputData[9];
	gam2=inputData[10];
	READ=inputData[11];
	wczytaj_sciany=inputData[12];
	granica_odczytu_scian=inputData[13];
	GAP=inputData[14];
}
VSH=round(1000*VSH)/1000;
VSH_K=round(1000*VSH_K)/1000;
P0=round(100*P0)/100;
P0_K=round(100*P0_K)/100;

cout<<VSH<<"  "<<VSH_K<<"  "<<P0<<"  "<<P0_K<<endl;




//SCIEZKI ZAPISU----------------------------------------------------------------

ofstream parametry;

ofstream R_X;
ofstream R_Y;
ofstream R_Z;

ofstream przebiegi;
ofstream v_prof;
ofstream T_prof;
ofstream H_prof;
ofstream ostatnie_wartosci;
ofstream hist_P_L;



parametry.open(nazwa_parametry);
if(!parametry)return 0;   

R_X.open(nazwa_RX);
if(!R_X)return 0; 

R_Y.open(nazwa_RY);
if(!R_X)return 0; 

R_Z.open(nazwa_RZ);
if(!R_Z)return 0; 

przebiegi.open(nazwa_przebiegi);
if(!przebiegi)return 0;

v_prof.open(nazwa_przekroje_v);
if(!v_prof)return 0;

T_prof.open(nazwa_przekroje_t);
if(!T_prof)return 0;

H_prof.open(nazwa_przekroje_h);
if(!H_prof)return 0;

ostatnie_wartosci.open(nazwa_ostatnie_wartosci);
if(!ostatnie_wartosci)return 0;

hist_P_L.open(nazwa_hist_P_L);
if(!hist_P_L)return 0;

//------------------------------------------------------------------------------
//READ DATA FROM FILE
//------------------------------------------------------------------------------
if(READ==1)
{
fstream plik_poloz1("input.txt", ios::in | ios::out );

TOTPX=0; TOTPY=0; TOTPZ=0;
        
double pomoc1, pomoc2, pomoc3, pomoc4, pomoc5, pomoc6, pomoc7, pomoc8, pomoc9;
	for (int i=0; i<NATOMS; i++) 
	{
	plik_poloz1 >> pomoc1 >> pomoc2 >> pomoc3 >> pomoc4 >> pomoc5 >> pomoc6 >> pomoc7 >> pomoc8 >> pomoc9;
	RX[i]=pomoc1;
	RY[i]=pomoc2;
	RZ[i]=pomoc3;
    VX[i]=pomoc4;
	VY[i]=pomoc5;
	VZ[i]=pomoc6;
	RX_eq[i]=pomoc1;
    RY_eq[i]=pomoc2;
    RZ_eq[i]=pomoc3;
	}
		
	for(int i=0;i<NATOMS;i++)
    {        
    TOTPX=TOTPX+VX[i];
    TOTPY=TOTPY+VY[i];
    TOTPZ=TOTPZ+VZ[i];
	}
	plik_poloz1.close();

if(ZerowaniePedu==1)
{    
	for (int i=0; i<NATOMS; i++)
    {  
    VX[i]=VX[i]-dziel(TOTPX,NATOMS);
    VY[i]=VY[i]-dziel(TOTPY,NATOMS);
    VZ[i]=VZ[i]-dziel(TOTPZ,NATOMS);   
    }
}

EKC=0;
for(int i=0;i<NATOMS;i++)
{
EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);
}
T=2*pow(3*NATOMS,-1)*EKC;
scale=T*pow(T0,-1);

SXP=0; SYP=0; SZP=0;
		for (int i=0; i<NATOMS; i++) 
        {    		
    		SXP=SXP+RX[i];
    		SYP=SYP+RY[i];
    		SZP=SZP+RZ[i];
        }
SXP=SXP*pow(NATOMS,-1); SYP=SYP*pow(NATOMS,-1); SZP=SZP*pow(NATOMS,-1);     

if(ZerowanieSrodkaM==1)
{
for (int i=0; i<NATOMS; i++) 
        {
            RX[i]=RX[i]-SXP;
            RY[i]=RY[i]-SYP;
            RZ[i]=RZ[i]-SZP;
            
            RX_eq[i]=RX[i];
    		RY_eq[i]=RY[i];
  		    RZ_eq[i]=RZ[i];
   		
       }
}
}
//------------------------------------------------------------------------------
if(READ==2)
{	
DBOX=BX;	

DSTX=BX*pow(NUNITX,-1);
DSTY=BY*pow(NUNITY,-1);
DSTZ=BZ*pow(NUNITZ,-1);

cax=DSTX;
cay=DSTY;
caz=DSTZ;

RX[0]=0;				RY[0]=0;					RZ[0]=0;
RX[1]=0;				RY[1]=DSTY*pow(2,-1);		RZ[1]=DSTZ*pow(2,-1);
RX[2]=DSTX*pow(2,-1);	RY[2]=0;					RZ[2]=DSTZ*pow(2,-1);
RX[3]=DSTX*pow(2,-1);	RY[3]=DSTY*pow(2,-1);		RZ[3]=0;

M=0;
GRAN=0;
        for (int i=0; i<NUNITZ; i++) //z
        {
            for (int j=0; j<NUNITY; j++) //y
            {
                for (int k=0; k<NUNITX; k++) //x
                {
                    for (int ij=0; ij<4; ij++) 
                    {
						RandA=rand(); if(RandA==0) {RandA=rand();} 
						RandB=rand(); if(RandB==0){RandB=rand();}
						
                    	GRAN=cos(2*M_PI*(RandA*pow(RAND_MAX,-1)))*pow(-2*log((RandB*pow(RAND_MAX,-1))),(1/2));
                        RX[ij+M]=RX[ij]+DSTX*(k)+0.05*GRAN;
                        
						RandA=rand(); if(RandA==0){RandA=rand();}
						RandB=rand(); if(RandB==0){RandB=rand();}
						
						GRAN=cos(2*M_PI*(RandA*pow(RAND_MAX,-1)))*pow(-2*log((RandB*pow(RAND_MAX,-1))),(1/2));
                        RY[ij+M]=RY[ij]+DSTY*(j)+0.05*GRAN;
                        
						RandA=rand(); if(RandA==0){RandA=rand();}
						RandB=rand(); if(RandB==0){RandB=rand();}
                    	
                        GRAN=cos(2*M_PI*(RandA*pow(RAND_MAX,-1)))*pow(-2*log((RandB*pow(RAND_MAX,-1))),(1/2));
                        RZ[ij+M]=RZ[ij]+DSTZ*(i)+0.05*GRAN;
                    }
                    M=M+4;
                }
            }
        }

		SXP=0; SYP=0; SZP=0;
		for (int i=0; i<NATOMS; i++) 
        {
            RX[i]=RX[i]+DSTX*pow(4,-1)-BX*pow(2,-1);
            RY[i]=RY[i]+DSTY*pow(4,-1)-BY*pow(2,-1);
            RZ[i]=RZ[i]+DSTZ*pow(4,-1)-BZ*pow(2,-1);
    		
    		SXP=SXP+RX[i];
    		SYP=SYP+RY[i];
    		SZP=SZP+RZ[i];
        }
        
  
SXP=SXP*pow(NATOMS,-1); SYP=SYP*pow(NATOMS,-1); SZP=SZP*pow(NATOMS,-1);     
for (int i=0; i<NATOMS; i++) 
        {
            RX[i]=RX[i]-SXP;
            RY[i]=RY[i]-SYP;
            RZ[i]=RZ[i]-SZP;
            
            RX_eq[i]=RX[i];
    		RY_eq[i]=RY[i];
    		RZ_eq[i]=RZ[i];
    		
        }
  SXP=0; SYP=0; SZP=0;      
        for (int i=0; i<NATOMS; i++) 
        {
           	SXP=SXP+RX[i];
    		SYP=SYP+RY[i];
    		SZP=SZP+RZ[i];    		
        }


//WCZYTYWANIE SCIAN-------------------------------------------------------------

if(wczytaj_sciany==1)
{
	fstream plik_poloz1("input.txt", ios::in | ios::out );
        
	double pomoc1, pomoc2, pomoc3, pomoc4, pomoc5, pomoc6, pomoc7, pomoc8, pomoc9;
	for (int i=0; i<NATOMS; i++) 
	{	
	plik_poloz1 >> pomoc1 >> pomoc2 >> pomoc3 >> pomoc4 >> pomoc5 >> pomoc6 >> pomoc7 >> pomoc8 >> pomoc9;
	if(i<granica_odczytu_scian)
	{
	RX[i]=pomoc1;
	RY[i]=pomoc2;
	RZ[i]=pomoc3;
    VX[i]=pomoc4;
	VY[i]=pomoc5;
	VZ[i]=pomoc6;
	RX_eq[i]=RX[i];
    RY_eq[i]=RY[i];
    RZ_eq[i]=RZ[i];
	}
	}

	SXP=0; SYP=0; SZP=0;      
        for (int i=0; i<NATOMS; i++) 
        {
           	SXP=SXP+RX[i];
    		SYP=SYP+RY[i];
    		SZP=SZP+RZ[i];
    		
        }
   

SXP=SXP*pow(NATOMS,-1); SYP=SYP*pow(NATOMS,-1); SZP=SZP*pow(NATOMS,-1);     

for (int i=0; i<NATOMS; i++) 
        {
            RX[i]=RX[i]-SXP;
            RY[i]=RY[i]-SYP;
            RZ[i]=RZ[i]-SZP;
            
            RX_eq[i]=RX[i];
    		RY_eq[i]=RY[i];
    		RZ_eq[i]=RZ[i];
    		
        }
  SXP=0; SYP=0; SZP=0;      
        for (int i=0; i<NATOMS; i++) 
        {
           	SXP=SXP+RX[i];
    		SYP=SYP+RY[i];
    		SZP=SZP+RZ[i];
    		
        }
	
}            


//------------------------------------------------------------------------------
        for (int i=0; i<wskA; i++) 
        {
            RZ[i]=RZ[i]-GAP;
    		RZ_eq[i]=RZ[i];
        }
        
        for (int i=wskB; i<NATOMS; i++) 
        {
            RZ[i]=RZ[i]+GAP;
    		RZ_eq[i]=RZ[i];
        }    

for (int i=0; i<0.5*NATOMS; i++) 
        {
            RZ[i]=RZ[i]-CENTGAP;
    		RZ_eq[i]=RZ[i];
        }
        
        for (int i=0.5*NATOMS; i<NATOMS; i++) 
        {
            RZ[i]=RZ[i]+CENTGAP;
    		RZ_eq[i]=RZ[i];
        }    
//------------------------------------------------------------------------------        


    for(int i=0;i<NATOMS;i++)
    {        
    
    VX[i]=0; VY[i]=0; VZ[i]=0;
       
    RandA=rand(); if(RandA==0) {RandA=rand();}
	RandB=rand(); if(RandB==0){RandB=rand();}
	   
    GRAN=cos(2*M_PI*(RandA*pow(RAND_MAX,-1)))*pow(-2*log((RandB*pow(RAND_MAX,-1))),(0.5));
    VX[i]=GRAN*T0;
    
    RandA=rand(); if(RandA==0){RandA=rand();}
	RandB=rand(); if(RandB==0) {RandB=rand();}
    
    GRAN=cos(2*M_PI*(RandA*pow(RAND_MAX,-1)))*pow(-2*log((RandB*pow(RAND_MAX,-1))),(0.5));
    VY[i]=GRAN*T0;
    
    RandA=rand(); if(RandA==0) {RandA=rand();}
	RandB=rand(); if(RandB==0) {RandB=rand();}
    
    GRAN=cos(2*M_PI*(RandA*pow(RAND_MAX,-1)))*pow(-2*log((RandB*pow(RAND_MAX,-1))),(0.5));
	VZ[i]=GRAN*T0;
    }

    
//------------------------------------------------------------------------------     
NATIN=wskA;
TOTPX=0; TOTPY=0; TOTPZ=0;        
    
	for (int i=0; i<wskA; i++)
    {  
    TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];  
    }
	  
    for (int i=0; i<wskA; i++)
    {  
    VX[i]=VX[i]-dziel(TOTPX,NATIN); 
	VY[i]=VY[i]-dziel(TOTPY,NATIN); 
	VZ[i]=VZ[i]-dziel(TOTPZ,NATIN);   
    }

TOTPX=0; TOTPY=0; TOTPZ=0;
for (int i=0; i<wskA; i++)
    {  	
	TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];
    } 

EKC=0;
for(int i=0;i<wskA;i++)
{
EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);
}
T=2*pow(3*NATIN,-1)*EKC;
scale=T*pow(T0,-1);

TOTPX=0; TOTPY=0; TOTPZ=0;
for (int i=0; i<wskA; i++)
    {  
    VX[i]=VX[i]*pow(scale,-0.5); 
	VY[i]=VY[i]*pow(scale,-0.5); 
	VZ[i]=VZ[i]*pow(scale,-0.5);   
	TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];
    }
//------------------------------------------------------------------------------
NATIN=NATOMS-2*wskA;
TOTPX=0; TOTPY=0; TOTPZ=0;        
    
	for (int i=wskA; i<wskB; i++)
    {  
    TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];  
    }
	  
    for (int i=wskA; i<wskB; i++)
    {  
    VX[i]=VX[i]-dziel(TOTPX,NATIN); 
	VY[i]=VY[i]-dziel(TOTPY,NATIN); 
	VZ[i]=VZ[i]-dziel(TOTPZ,NATIN);   
    }

TOTPX=0; TOTPY=0; TOTPZ=0;
for (int i=wskA; i<wskB; i++)
    {  	
	TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];
    }

EKC=0;
for(int i=wskA;i<wskB;i++)
{
EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);
}
T=2*pow(3*NATIN,-1)*EKC;
scale=T*pow(T0,-1);

TOTPX=0; TOTPY=0; TOTPZ=0;
for (int i=wskA; i<wskB; i++)
    {  
    VX[i]=VX[i]*pow(scale,-0.5); 
	VY[i]=VY[i]*pow(scale,-0.5); 
	VZ[i]=VZ[i]*pow(scale,-0.5);   

	TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];
    }
//------------------------------------------------------------------------------
NATIN=wskA;
TOTPX=0; TOTPY=0; TOTPZ=0;        
    
	for (int i=wskB; i<NATOMS; i++)
    {  
    TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];  
    }
	  
    for (int i=wskB; i<NATOMS; i++)
    {  
    VX[i]=VX[i]-dziel(TOTPX,NATIN); 
	VY[i]=VY[i]-dziel(TOTPY,NATIN); 
	VZ[i]=VZ[i]-dziel(TOTPZ,NATIN);   
    }


TOTPX=0; TOTPY=0; TOTPZ=0;
for (int i=wskB; i<NATOMS; i++)
    {  	
	TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];
    }


EKC=0;
for(int i=wskB;i<NATOMS;i++)
{
EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);
}
T=2*pow(3*NATIN,-1)*EKC;
scale=T*pow(T0,-1);

TOTPX=0; TOTPY=0; TOTPZ=0;
for (int i=wskB; i<NATOMS; i++)
    {  
    VX[i]=VX[i]*pow(scale,-0.5); 
	VY[i]=VY[i]*pow(scale,-0.5); 
	VZ[i]=VZ[i]*pow(scale,-0.5);  

	TOTPX=TOTPX+VX[i]; TOTPY=TOTPY+VY[i]; TOTPZ=TOTPZ+VZ[i];
    }
  
//------------------------------------------------------------------------------



EKC=0;
for(int i=0;i<NATOMS;i++)
{
EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);
}
T=2*pow(g,-1)*EKC;
scale=T*pow(T0,-1);



	if(StartVelocities==0)
	{
		for (int i=0; i<NATOMS; i++)
    	{  
    	VX[i]=0; 
		VY[i]=0; 
		VZ[i]=0;   
    	}
	}


}

//VOLUME------------------------------------------------------------------------

DIMZ=fabs(RZ[NATOMS-1]-RZ[0]);

PDIM1=0;
for (int gg=wskA-wskA*pow(2*LEY,-1); gg<wskA; gg++)
{
	PDIM1=PDIM1+RZ[gg];
}
PDIM1=PDIM1*pow(2*LEY,-1);

PDIM2=0;
for (int gg=wskB; gg<wskB+wskA*pow(2*LEY,-1); gg++)
{
	PDIM2=PDIM2+RZ[gg];
}
PDIM2=PDIM2*pow(2*LEY,-1);

DIMZ2=fabs(PDIM2-PDIM1);
DIMZ3=fabs(RZ[wskB]-RZ[wskA]);

Volum1=BX*BX*DIMZ;
Volum2=BX*BX*DIMZ2;
Volum4=BX*BX*DIMZ3;

EKZ=0;
for(int i=0;i<NATOMS;i++) {EKZ=EKZ+0.5*(VZ[i]*VZ[i]);}
EKZ=EKZ*pow(NATOMS,-1);
EKZ=EKZ*2*NATOMS*pow(Volume,-1);
//------------------------------------------------------------------------------
//STATE CALCULATION
//------------------------------------------------------------------------------
counter=0;
                    
        for(int i=0;i<NATOMS;i++)
		{
            NLL=-1;
            marker[i]=0;
            for(int j=0;j<NATOMS;j++)
			{
                            
            X=0; Y=0; Z=0; R2=0; R=0; RI=0;
			X=RX[i]-RX[j]; Y=RY[i]-RY[j]; Z=RZ[i]-RZ[j];
    
			X=X-XPERIOD*(DBOX*round(X*pow(DBOX,-1)));
			Y=Y-YPERIOD*(DBOX*round(Y*pow(DBOX,-1)));
			Z=Z-ZPERIOD*(DBOZ*round(Z*pow(DBOZ,-1)));
    
			R2=X*X+Y*Y+Z*Z;
                            
                if ((R2<((RCUT+SKIN)*(RCUT+SKIN))) && (i!=j)) 
				{
                    NLL++;
                    list[i][NLL]=j;
                }
            }
            marker[i]=NLL+1;     
			
        }
//counter=counter+1;        



PXX=0; PYY=0; PZZ=0; PXY=0, PYX=0; PXZ=0; PZX=0; PZY=0; PYZ=0;
PXX_mid=0; PYY_mid=0; PZZ_mid=0;
FUP2=0; FDOWN2=0;
for (int i=0;i<NATOMS;i++)
{
	
FX[i]=0.0;
FY[i]=0.0;
FZ[i]=0.0;
U[i]=0.0;
ULJ[i]=0.0;
UPH[i]=0.0;

EUP=0.0;
FXP=0.0;
FYP=0.0;
FZP=0.0;  

	for(int j=0;j<NATOMS;j++)
	{
		        
        X=0; Y=0; Z=0; R2=0; R=0; RI=0;
		X=RX[i]-RX[j]; Y=RY[i]-RY[j]; Z=RZ[i]-RZ[j];
    
		X=X-XPERIOD*(DBOX*round(X*pow(DBOX,-1)));
		Y=Y-YPERIOD*(DBOX*round(Y*pow(DBOX,-1)));
		Z=Z-ZPERIOD*(DBOZ*round(Z*pow(DBOZ,-1)));

		R2=X*X+Y*Y+Z*Z;
		R=pow(R2,0.5);
		RI=pow(R,-1);
	
		if(i!=j)
		{
		//----------------------------------------------------------------------	
			
			if(R2<RCUT*RCUT)
			{
			//------------------------------------------------------------------	
				if(DIFINT==1)
				{
				if((i<wskA) && (j<wskA)) 								{a=a1; b=1;}	
				else if((i>=wskB) && (j>=wskB)) 						{a=a1; b=1;}
				else if((i>=wskA) && (j>=wskA) && (i<wskB) && (j<wskB)) {a=1; b=1;}
				else 													{a=a2; b=1;}
				}
				else {a=1; b=1;}
				
			//------------------------------------------------------------------				
			EUP=4*a*(pow(RI,N1)-b*pow(RI,N2));
			EUP=EUP+shift;	  
      		
			FXP=4*a*(N1*X*pow(RI,N1+2)-b*N2*X*pow(RI,N2+2));
			FYP=4*a*(N1*Y*pow(RI,N1+2)-b*N2*Y*pow(RI,N2+2));
			FZP=4*a*(N1*Z*pow(RI,N1+2)-b*N2*Z*pow(RI,N2+2));
			
				if((R2<R3*R3 && i<wskA && j<wskA) || (R2<R3*R3 && i>=wskB && j>=wskB))
				{
				FXP=FXP-PAIR_TET*WR3*KR3*pow(R,WR3-2)*X;
				FYP=FYP-PAIR_TET*WR3*KR3*pow(R,WR3-2)*Y;
				FZP=FZP-PAIR_TET*WR3*KR3*pow(R,WR3-2)*Z;
				
				EUP=EUP+PAIR_TET*KR3*pow(R,WR3);
				}
					
			}   
			else if(R2>=RCUT*RCUT) {EUP=0.0; FXP=0.0; FYP=0.0; FZP=0.0;}   

			FX[i]+=FXP; FY[i]+=FYP; FZ[i]+=FZP; U[i]+=0.5*EUP; 
		//----------------------------------------------------------------------
		}
    }
	ULJ[i]=U[i];   
}

EPC_walls=0;
for(int i=0;i<wskA;i++) {EPC_walls=EPC_walls+U[i];}
for(int i=wskB;i<NATOMS;i++) {EPC_walls=EPC_walls+U[i];}

EPC_w_down=0;
for(int i=0;i<wskA;i++) {EPC_w_down=EPC_w_down+U[i];}

EPC_w_up=0;
for(int i=wskB;i<NATOMS;i++) {EPC_w_up=EPC_w_up+U[i];}

EPC_inter=0;
for(int i=wskA;i<wskB;i++) {EPC_inter=EPC_inter+U[i];}

EPC=0;
for(int i=0;i<NATOMS;i++) {EPC=EPC+U[i];}

EKC=0;
g_in=0;
probny_counter=0;
for(int i=0;i<wskA;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);probny_counter=probny_counter+1;}
g_in=3*(wskA);
T_down=2*pow(g_in,-1)*EKC;

EKC=0;
g_in=0;
probny_counter=0;
for(int i=wskB;i<NATOMS;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);probny_counter=probny_counter+1;}
g_in=3*(wskA);
T_up=2*pow(g_in,-1)*EKC;

EKC=0;
g_in=0;
for(int i=wskA;i<wskB;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
g_in=3*(NATOMS-2*wskA);
T_mid=2*pow(g_in,-1)*EKC;
cout<<EKC<<endl;

EKC=0;
g_in=0;
for(int i=wskA;i<wskB;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
g_in=3*(NATOMS-2*wskA);
T_mid_2=2*pow(g_in,-1)*EKC;

EKC=0;
g_in=0;
for(int i=0;i<wskA;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
for(int i=wskB;i<NATOMS;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
g_in=3*(2*wskA);
T_walls=2*pow(g_in,-1)*EKC;
cout<<EKC<<endl;

EKC=0;
for(int i=0;i<NATOMS;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
T=2*pow(g,-1)*EKC;

EKC=0;
for(int i=0;i<NATOMS;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}

EKC_walls=0;
for(int i=0;i<wskA;i++) {EKC_walls=EKC_walls+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
for(int i=wskB;i<NATOMS;i++) {EKC_walls=EKC_walls+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}

EKC_w_down=0;
for(int i=0;i<wskA;i++) {EKC_w_down=EKC_w_down+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}

EKC_w_up=0;
for(int i=wskB;i<NATOMS;i++) {EKC_w_up=EKC_w_up+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}

EKC_inter=0;
for(int i=wskA;i<wskB;i++) {EKC_inter=EKC_inter+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}

calka=(EPC+EKC)*pow(NATOMS,-1);;   
calkaBEGIN=calka;
//------------------------------------------------------------------------------
//HALF STEP
//------------------------------------------------------------------------------
for(int i=0;i<NATOMS;i++)
{
RXR[i]=RX[i];
RYR[i]=RY[i];
RZR[i]=RZ[i];	
	
VXH1[i]=VX[i]+0.5*DT*FX[i]; VYH1[i]=VY[i]+0.5*DT*FY[i]; VZH1[i]=VZ[i]+0.5*DT*FZ[i];
VXH2[i]=VXH1[i]; VYH2[i]=VYH1[i]; VZH2[i]=VZH1[i];
}

//------------------------------------------------------------------------------

CPX=0;
CPY=0;
CPZ=0;

CPXA=0;
CPYA=0;
CPZA=0;

CPXB=0;
CPYB=0;
CPZB=0;

CPXC=0;
CPYC=0;
CPZC=0;

SMX=0;
SMY=0;
SMZ=0;

SMXA=0;
SMYA=0;
SMZA=0;

SMXB=0;
SMYB=0;
SMZB=0;

SMXC=0;
SMYC=0;
SMZC=0;

CMPX=0;
CMPY=0;
CMPZ=0;

for (int ii=0;ii<NATOMS;ii++)
{
	if(ii<wskA)
	{
	CPXA=CPXA+VX[ii];
	CPYA=CPYA+VY[ii];
	CPZA=CPZA+VZ[ii];	
	
	SMXA=SMXA+RX[ii];
	SMYA=SMYA+RY[ii];
	SMZA=SMZA+RZ[ii];
	}
	else if(ii>=wskB)
	{
	CPXC=CPXC+VX[ii];
	CPYC=CPYC+VY[ii];
	CPZC=CPZC+VZ[ii];	
	
	SMXC=SMXC+RX[ii];
	SMYC=SMYC+RY[ii];
	SMZC=SMZC+RZ[ii];
	}
	else
	{
	CPXB=CPXB+VX[ii];
	CPYB=CPYB+VY[ii];
	CPZB=CPZB+VZ[ii];
	
	SMXB=SMXB+RX[ii];
	SMYB=SMYB+RY[ii];
	SMZB=SMZB+RZ[ii];	
	}
	
CPX=CPX+VX[ii];
CPY=CPY+VY[ii];
CPZ=CPZ+VZ[ii];

SMX=SMX+RX[ii];
SMY=SMY+RY[ii];
SMZ=SMZ+RZ[ii];

CMPX=CMPX+(VZ[ii]*RY[ii]-VY[ii]*RZ[ii]);
CMPY=CMPY+(VX[ii]*RZ[ii]-VZ[ii]*RX[ii]);
CMPZ=CMPZ+(VY[ii]*RX[ii]-VX[ii]*RY[ii]);
}

SMX=SMX*pow(NATOMS,-1);
SMY=SMY*pow(NATOMS,-1);
SMZ=SMZ*pow(NATOMS,-1);

SMXA=SMXA*pow(NATOMS,-1);
SMYA=SMYA*pow(NATOMS,-1);
SMZA=SMZA*pow(NATOMS,-1);

SMXB=SMXB*pow(NATOMS,-1);
SMYB=SMYB*pow(NATOMS,-1);
SMZB=SMZB*pow(NATOMS,-1);

SMXC=SMXC*pow(NATOMS,-1);
SMYC=SMYC*pow(NATOMS,-1);
SMZC=SMZC*pow(NATOMS,-1);
//------------------------------------------------------------------------------
//ZAPIS PARAMETROW SYMULACJI
//------------------------------------------------------------------------------

parametry.precision(16);

//number and index of particle-------------------------------------------------- 

parametry<</*1*/NATOMS<<endl;
parametry<</*2*/wskA<<endl;
parametry<</*3*/wskB<<endl;

//run parameters----------------------------------------------------------------

parametry<</*4*/RUN_TIME<<endl;
parametry<</*5*/RUN<<endl;
parametry<</*6*/DT<<endl;
parametry<</*7*/cleaning<<endl;
parametry<</*8*/JUMP<<endl;
parametry<</*9*/JUMP2<<endl;

//thermodynamics----------------------------------------------------------------

parametry<</*10*/mass<<endl;
parametry<</*11*/KB<<endl;
parametry<</*12*/T0<<endl;
parametry<</*13*/T0_up<<endl;
parametry<</*14*/T0_down<<endl;
parametry<</*15*/g<<endl;
parametry<</*16*/rho<<endl;
parametry<</*17*/Volume<<endl;
parametry<</*18*/P0<<endl;
parametry<</*19*/Qp<<endl;
parametry<</*20*/tauH<<endl;
parametry<</*21*/tau<<endl;
parametry<</*22*/Q_in<<endl;
parametry<</*23*/Q<<endl;
parametry<</*24*/RCUT<<endl;
parametry<</*25*/shift<<endl;
parametry<</*26*/SKIN<<endl;

//geometry----------------------------------------------------------------------

parametry<</*27*/DBOX<<endl;
parametry<</*28*/BX<<endl;
parametry<</*29*/BZ<<endl;
parametry<</*30*/GAP<<endl;
parametry<</*31*/CENTGAP<<endl;

//logic variables---------------------------------------------------------------

parametry<</*32*/READ<<endl;
parametry<</*33*/VSH<<endl;
parametry<</*34*/SHEAR<<endl;
parametry<</*35*/BAR_CHOICE<<endl;
parametry<</*36*/StartVelocities<<endl;
parametry<</*37*/XPERIOD<<endl;
parametry<</*38*/YPERIOD<<endl;
parametry<</*39*/ZPERIOD<<endl;
parametry<</*40*/TERM<<endl;
parametry<</*41*/termWall<<endl;
parametry<</*42*/TETF<<endl;
parametry<</*43*/BAR<<endl;
parametry<</*44*/CONSTR<<endl;
parametry<</*45*/DIFINT<<endl;
parametry<</*46*/DIFCON<<endl;
parametry<</*47*/termX<<endl;
parametry<</*48*/termY<<endl;
parametry<</*49*/termZ<<endl;
parametry<</*50*/cax<<endl;
parametry<</*51*/cay<<endl;
parametry<</*52*/caz<<endl;
parametry<</*53*/omega<<endl;
parametry<</*54*/beta<<endl;
parametry<</*55*/gam1<<endl;
parametry<</*56*/gam2<<endl;
parametry<</*57*/P0<<endl;
parametry<</*58*/P0_K<<endl;
parametry<</*59*/VSH<<endl;
parametry<</*60*/VSH_K<<endl;

parametry.flush();
parametry.close();

//------------------------------------------------------------------------------
deltay=30*pow(PGAP2,-1);

	for(int kk=0;kk<PGAP2;kk++)	
	{
	POZSURF[kk]=deltay*(kk-0.5*PGAP2);
	kinet[kk]=0;
	timkinet[kk]=0;
	}
	
//------------------------------------------------------------------------------




cout<<CPXC<<"   "<<CPYC<<"   "<<CPZC<<endl;
cout<<"*****************************************************************"<<endl;
cout<<" "<<endl;

cout<<"   PROGRRAM - N-DEP"<<endl;

cout<<" "<<endl;

cout<<"   N:           "<<NATOMS<<endl;
cout<<"   wskA:        "<<wskA<<endl;
cout<<"   wskB:        "<<wskB<<endl;
cout<<"   RUN:         "<<RUN<<endl;
cout<<"   DT:          "<<DT<<endl;
cout<<"   READ:        "<<READ<<endl;

cout<<" "<<endl;
cout<<"-----------------------------------------------------------------"<<endl;
cout<<" "<<endl;

cout<<"   rho:         "<<rho<<endl;
cout<<"   T0:          "<<T0<<endl;
cout<<"   P0:          "<<P0<<endl;
cout<<"   Vsh:		   "<<VSH<<endl;

cout<<" "<<endl;
cout<<"-----------------------------------------------------------------"<<endl;
cout<<" "<<endl;

cout<<"   RCUT:        "<<RCUT<<endl;
cout<<"   shift:       "<<shift<<endl;
cout<<"   skin:        "<<SKIN<<endl;

cout<<" "<<endl;
cout<<"-----------------------------------------------------------------"<<endl;
cout<<" "<<endl;

cout<<"   GAP:         "<<GAP<<endl;
cout<<"   StartVel:    "<<StartVelocities<<endl;
cout<<"   XPERIOD:     "<<XPERIOD<<endl;
cout<<"   YPERIOD:     "<<YPERIOD<<endl;
cout<<"   ZPERIOD:     "<<ZPERIOD<<endl;
cout<<"   S:		   "<<BX*BY<<endl;

cout<<" "<<endl;
cout<<"-----------------------------------------------------------------"<<endl;
cout<<" "<<endl;

cout<<"   cleaning:    "<<cleaning<<endl;
cout<<"   JUMP:        "<<JUMP<<endl;

cout<<" "<<endl;	
cout<<"   VSH: "<<VSH<<" VSH_K: "<<VSH_K<<" P0: "<<P0<<" P0_K: "<<P0_K<<" RUN: "<<RUN<<" cl: "<<cleaning<<" a1: "<<a1<<" b4: "<<BARO4<<" b5: "<<BARO5<<" READ: "<<READ<<endl;

cout<<" "<<endl;
cout<<"*****************************************************************"<<endl;

//------------------------------------------------------------------------------
//MAIN PROGRAMME
//------------------------------------------------------------------------------

t1=clock();

time_t rawtime;	
struct tm * timeinfo;
time (&rawtime);
timeinfo = localtime (&rawtime);
cout<<asctime(timeinfo)<<endl;

for (int n=0;n<RUN;n++)
{	


if(n==t_test)
{
	t2=clock();
	t3=(RUN/t_test)*(t2-t1);
	
	
	cout<<"czas rozpoczecia petli glownej [s]: "<<t1/sekunda<<endl;
	cout<<"czas obliczania pierszej serii [s]: "<<(t2-t1)/sekunda<<endl;
	cout<<"Szacowany czas trwania obliczen [h]: "<<t3/godzina/1000<<endl;
	
}


fstream plik_time(nazwa_time, ios::out );
    if (plik_time.good()) 
	{  
	plik_time<<asctime(timeinfo)<<endl;
	plik_time<<"czas rozpoczecia petli glownej [s]: "<<t1/sekunda<<endl;
	plik_time<<"czas obliczania pierszej serii [s]: "<<(t2-t1)/sekunda<<endl;
	plik_time<<"Szacowany czas trwania obliczen [h]: "<<t3/godzina/1000<<endl;
    //zamkniecie pliku
    plik_time.close();
    }


if(n==0.5*cleaning)
{
	
P0=P0_K;
VSH=VSH_K; 
automat_sliding=0;
automat_pressing=0;
deltaV=0;
deltaP=0;
pressing=0;
sliding=0;
	
}



if(n==cleaning)
{	
GET=1;	
kincount=0;
CNT=0;
P0=P0_K;
VSH=VSH_K; //poczatkowa predkosc scinania
//scale_vel_w=0;
Qp=Qp;	
	
DT=DT;	
tim=0;	
DENOM=1;
liczsumk=0;

sumVW=0;

sumNEWP1=0;
sumNEWP2=0;

sumP=0; 
SUMFRIC1=0;
SUMFRIC2=0; 
SUMFRIC3=0;
SUMFRIC4=0;

SUMFRIC5=0;
SUMFRIC6=0;
SUM_EP=0;
SUM_T=0;

SUM_T_up=0;
SUM_T_down=0;
SUM_T_mid=0;
SUM_T_mid_2=0;
SUM_T_walls=0;
SUM_T_pec=0;

SUM_P=0;
sumPXX=0;
sumPYY=0;
sumPZZ=0;
sumPTX=0;
sumPTY=0;
sumPTZ=0;
sumPVZ=0;
sumDIM=0;	

sumPYZ=0;
sumPZY=0;
sumPXZ=0;
sumPZX=0;
sumPYX=0;
sumPXY=0;

sumPXX_mid=0;
sumPYY_mid=0;
sumPZZ_mid=0;
sumPYZ_mid=0;

sumPT1=0; 
sumPT2=0;

sumPWX=0;
sumPWY=0;

SUM_PXX=0;
SUM_PYY=0;
SUM_PZZ=0;
SUM_PXY=0; 
SUM_PYX=0;
SUM_PZX=0; 
SUM_PXZ=0; 
SUM_PYZ=0;
SUM_PZY=0;
SUM_FDOWN2=0;
SUM_FUP2=0;
SUM_FBXD=0;
SUM_FBXU=0;
SUM_FBYD=0;
SUM_FBYU=0;

sumFTET1=0; 
sumFTET2=0;

sumEPC_inter=0;
sumEPC_w_down=0;
sumEPC_w_up=0;
meanEPC_inter=0; 
meanEPC_w_down=0;
meanEPC_w_up=0;

SUM_DR=0;
meanDR=0;


sumDEW1=0;
sumDEW2=0;

sumSMX=0;
sumSMY=0;
sumSMZ=0;
sumSMXA=0;
sumSMYA=0; 
sumSMZA=0; 
sumSMXB=0;
sumSMYB=0; 
sumSMZB=0; 
sumSMXC=0; 
sumSMYC=0; 
sumSMZC=0;

sumL=0;

sumMI=0;
sumMI1=0;
sumMI2=0;
sumMI3=0;
sumMI4=0;
sumMI5=0;
sumMI6=0;

VWU=0;
VWD=0;

sumTX=0;
sumTY=0;
sumTZ=0;

sumVolum1=0;
sumVolum2=0;
sumVolum3=0;

sumPVIR1=0;
sumPVIR2=0;
sumPVIR3=0;
sumPVIR4=0;
sumSPVIR1=0;
sumSPVIR2=0;
sumSPVIR3=0;

sumVAP=0;
meanVAP=0; 

sumVAP2=0;
meanVAP2=0;
sumVAP3=0;
meanVAP3=0;
sumVAP4=0;
meanVAP4=0;
sumVAP5=0;
meanVAP5=0;
sumVAP10=0;
meanVAP10=0;
sumVAP20=0;
meanVAP20=0;

sumPVIR_mid1=0;
sumPVIR_mid2=0;
sumPVIR_mid3=0;

meanPVIR_mid1=0;
meanPVIR_mid2=0;
meanPVIR_mid3=0;

sumPW=0;

distA=0; 
distC=0;

sumShFA=0;
sumShFC=0;

for(int j=0;j<TGAP;j++)
{
SUM_RTZ[j]=0;
SUM_RTZ_2[j]=0;

SUM_RTZ_X[j]=0;
SUM_RTZ_Y[j]=0;
SUM_RTZ_Z[j]=0;

SUM_RTVX[j]=0;
SUM_RTVY[j]=0;
SUM_RTVZ[j]=0;

tabDenomRTZ[j]=0;	
tabDenomRTZ_2[j]=0; 
tabDenomRTZ_X[j]=0; 	
tabDenomRTZ_Y[j]=0; 	
tabDenomRTZ_Z[j]=0;
tabDenomRVX[j]=0; 	
tabDenomRVY[j]=0;	
tabDenomRVY_pec[j]=0; 	
tabDenomRVZ[j]=0; 

histVXA[j]=0;
histVYA[j]=0; 
histVZA[j]=0; 
histVA[j]=0;
histVXB[j]=0; 
histVYB[j]=0; 
histVZB[j]=0; 
histVB[j]=0;
histVXC[j]=0; 
histVYC[j]=0; 
histVZC[j]=0; 
histVC[j]=0;

}



}

//------------------------------------------------------------------------------
if(zmiana_warunkow_T==1)
{	
	if(n==n_0)
	{
		scale_vel_w=2;
	}
}
//------------------------------------------------------------------------------
if(automat_sliding==1)
{
	if(VSH<VSH_K)
	{
		sliding=1;
	}
	else if(VSH>VSH_K)
	{
		sliding=2;
	}
}

if(sliding==1 && n>n_0 && n%n_0==0 && n<0.5*cleaning && VSH<VSH_K)
{
	
	VSH=VSH+double(deltaV*int(100*(rand()*pow(RAND_MAX,-1))+0.5))/double(10000);
	
}
if(sliding==2 && n>n_0 && n%n_0==0 && n<0.5*cleaning && VSH>VSH_K)
{
	
	VSH=VSH-double(deltaV*int(100*(rand()*pow(RAND_MAX,-1))+0.5))/double(10000);
	
}

if(VSH==VSH_K)
{
	
	if(flaga_AS==0)
	{
	automat_sliding=0;
	sliding=0;
	cout<<"Automatyczne rozpedzanie WYLACZONE!"<<endl;
	}
	
flaga_AS=1;
}
//------------------------------------------------------------------------------
if(automat_pressing==1)
{
	if(P0<P0_K)
	{
		pressing=1;
	}
	else if(P0>P0_K)
	{
		pressing=2;
	}
}



if(pressing==1 && n>n_0 && n%n_0==0 && n<0.5*cleaning && P0<P0_K)
{
	P0=P0+double(deltaP*int(100*(rand()*pow(RAND_MAX,-1))+0.5))/double(100);	
}
if(pressing==2 && n>n_0 && n%n_0==0 && n<0.5*cleaning && P0>P0_K)
{
	P0=P0-double(deltaP*int(100*(rand()/RAND_MAX)+0.5))/double(100);
}

if(P0==P0_K)
{
	
	if(flaga_AP==0)
	{
	automat_pressing=0;
	pressing=0;
	cout<<"Automatyczne sciskanie WYLACZONE!"<<endl;
	}
	
flaga_AP=1;
}

//------------------------------------------------------------------------------
if(jump_P0==1)
{
	if(n==6000) {P0=P0-deltaJumpP;  cleaning=n+1;}
	if(n==11000) {P0=P0+deltaJumpP; cleaning=n+1;}
}
//------------------------------------------------------------------------------
tim=tim+DT;

daik=1./(tim*BX*BX);
daic=1./(2*BX*BX*(tim/DT));

//obliczanie objetosci uk³adu---------------------------------------------------
POSZ1=0; POSZ2=0;
for (int i=0;i<wskA;i++)
{
	POSZ1=POSZ1+RZ[i]*pow(wskA,-1);
}
for (int i=wskB;i<NATOMS;i++)
{
	POSZ2=POSZ2+RZ[i]*pow(wskA,-1);
}

LENGTH=fabs(POSZ2-POSZ1);
sumL=sumL+LENGTH;

DIMZ=fabs(RZ[NATOMS-1]-RZ[0]);
	
PDIM1=0;
for (int gg=wskA-wskA*pow(2*LEY,-1); gg<wskA; gg++)
{
	PDIM1=PDIM1+RZ[gg];
}
PDIM1=PDIM1*pow(2*LEY,-1);

PDIM2=0;
for (int gg=wskB; gg<wskB+wskA*pow(2*LEY,-1); gg++)
{
	PDIM2=PDIM2+RZ[gg];
}
PDIM2=PDIM2*pow(2*LEY,-1);

DIMZ2=fabs(PDIM2-PDIM1);
DIMZ3=fabs(RZ[wskB]-RZ[wskA]);

Volum1=BX*BX*DIMZ;
Volum2=BX*BX*DIMZ2;
Volum3=BX*BX*LENGTH;
Volum4=BX*BX*DIMZ3;

sumVolum1=sumVolum1+Volum1;
sumVolum2=sumVolum2+Volum2;
sumVolum3=sumVolum3+Volum3;

meanVolum1=sumVolum1*pow(DENOM,-1);
meanVolum2=sumVolum2*pow(DENOM,-1);
meanVolum3=sumVolum3*pow(DENOM,-1);

//NEWTON EQUATION - VERTEL LEAP FROG--------------------------------------------	
//left right lists--------------------------------------------------------------
for(int i=0;i<NATOMS;i++)
{
	LX[NATOMS]=0;
	LY[NATOMS]=0;
	//1-mniejszy od 0; 2-wiekszy od 0;
	if (RX[i]<0)
	{
	LX[i]=1;	
	}
	else if(RX>0)
	{
	LX[i]=2;	
	}
	
	if (RY[i]<0)
	{
	LY[i]=1;	
	}
	else if(RY>0)
	{
	LY[i]=2;	
	}

}

//verlet list-------------------------------------------------------------------
int nmaksimum;
nmaksimum=NATOMS;
if(counter==1)
{
	//cout<<n<<endl;
counter=0;	
        for(int i=0;i<nmaksimum;i++)
		{
			DRX[i]=0;
			DRY[i]=0;
			DRZ[i]=0;
	
            NLL=-1;
            marker[i]=0;
            licznikomp=0;
            

            for(int jj=0;jj<nmaksimum;jj++)
			{
            licznikomp=licznikomp+1;                
            X=0; Y=0; Z=0; R2=0; R=0; RI=0;
			X=RX[i]-RX[jj]; Y=RY[i]-RY[jj]; Z=RZ[i]-RZ[jj];
    
			X=X-XPERIOD*(DBOX*round(X*pow(DBOX,-1)));
			Y=Y-YPERIOD*(DBOX*round(Y*pow(DBOX,-1)));
			Z=Z-ZPERIOD*(DBOZ*round(Z*pow(DBOZ,-1)));
    
			R2=X*X+Y*Y+Z*Z;
                            
                if ((R2<((RCUT+SKIN)*(RCUT+SKIN))) && (i!=jj)) 
				{
                    NLL++;
                    list[i][NLL]=jj;
                }
            }
            marker[i]=NLL+1;     
        }
               	
}
//------------------------------------------------------------------------------
PXX=0; PYY=0; PZZ=0; PXY=0; PYX=0; PXZ=0; PZX=0; PZY=0; PYZ=0; PVIR=0;	
PXX_mid=0; PYY_mid=0; PZZ_mid=0; PYZ_mid=0;
FUP2=0; FDOWN2=0; FXLEFT=0; FXRIGHT=0; FYLEFT=0; FYRIGHT=0;
F_up=0; F_down=0; FBXD=0; FBYD=0; FBXU=0; FBYU=0;

FP1U=0; FP1D=0; FP2U=0; FP2D=0; FP3U=0; FP3D=0; FP4U=0; FP4D=0; FP5U=0; FP5D=0; 
FP6U=0; FP6D=0; FP7U=0; FP7D=0; FP8U=0; FP8D=0; FP9U=0; FP9D=0; FP10U=0; FP10D=0;
FP11U=0; FP11D=0; FP12U=0; FP12D=0; FP13U=0; FP13D=0; FP14U=0; FP14D=0; FP15U=0; FP15D=0; 
FP16U=0; FP16D=0; FP17U=0; FP17D=0; FP18U=0; FP18D=0; FP19U=0; FP19D=0; FP20U=0; FP20D=0;
	

UAA=0;UA=0;UBB=0;UAB=0;UBC=0;UCC=0;UC=0;
FAAX=0;FAX=0;FBBX=0;FABX=0;FBCX=0;FCCX=0;FCX=0;
FAAY=0;FAY=0;FBBY=0;FABY=0;FBCY=0;FCCY=0;FCY=0;
FAAZ=0;FAZ=0;FBBZ=0;FABZ1=0;FABZ2=0;FBCZ1=0;FBCZ2=0;FCCZ=0;FCZ=0;
FTET1=0; FTET2=0;
ucount1=0;
ucount2=0;
count1=0; count2=0; count3=0; count4=0;
count5=0; count6=0; count7=0; count8=0;

int j;


		for (int i=0;i<NATOMS;i++)//single particle loop
		{	
		FX[i]=0.0;	FY[i]=0.0;	FZ[i]=0.0;	U[i]=0.0; licznikomp=0;

			for(int jj=0;jj<marker[i];jj++)//two-body loop
			{
			licznikomp=licznikomp+1;
			j=list[i][jj];
			if(i>j)
			{
        	
			EUP=0.0;	FXP=0.0;	FYP=0.0;	FZP=0.0;  
			X=0; Y=0; Z=0; R2=0; R=0; RI=0;
			X=RX[i]-RX[j]; Y=RY[i]-RY[j]; Z=RZ[i]-RZ[j];
    
			X=X-XPERIOD*(DBOX*round(X*pow(DBOX,-1)));
			Y=Y-YPERIOD*(DBOX*round(Y*pow(DBOX,-1)));
			Z=Z-ZPERIOD*(DBOZ*round(Z*pow(DBOZ,-1)));
    
			R2=X*X+Y*Y+Z*Z;
			R=pow(R2,0.5);
			RI=pow(R,-1);
			
		    //------------------------------------------------------------------	
			
			if(R2<RCUT*RCUT)
			{
			//------------------------------------------------------------------	
				if(DIFINT==1)
				{
				if((i<wskA) && (j<wskA)) 								{a=a1; b=1;} //A-A
				else if((i>=wskB) && (j>=wskB)) 						{a=a1; b=1;} //C-C
				else if((i>=wskA) && (j>=wskA) && (i<wskB) && (j<wskB)) {a=1; b=1;} //B-B
				else 													{a=a2; b=1;}//odpowiada oddzialywaniu sciana-plyn
				}
				else {a=1; b=1;}
				
			//------------------------------------------------------------------							
      		EUP=4*a*(pow(RI,N1)-b*pow(RI,N2));
			EUP=EUP+shift;
      		
			FXP=4*a*(N1*X*pow(RI,N1+2)-b*N2*X*pow(RI,N2+2));
			FYP=4*a*(N1*Y*pow(RI,N1+2)-b*N2*Y*pow(RI,N2+2));
			FZP=4*a*(N1*Z*pow(RI,N1+2)-b*N2*Z*pow(RI,N2+2));    
			
			if((R2<R3*R3 && i<wskA && j<wskA) || (R2<R3*R3 && i>=wskB && j>=wskB))
				{
				FXP=FXP-PAIR_TET*WR3*KR3*pow(R,WR3-2)*X;
				FYP=FYP-PAIR_TET*WR3*KR3*pow(R,WR3-2)*Y;
				FZP=FZP-PAIR_TET*WR3*KR3*pow(R,WR3-2)*Z;
				
				EUP=EUP+PAIR_TET*KR3*pow(R,WR3);
				}
		
			//skladowe tensora naprezen-----------------------------------------
			
			PXX+=2*FXP*X;
			PYY+=2*FYP*Y;
			PZZ+=2*FZP*Z;  
				if((i>=wskA) && (j>=wskA) && (i<wskB) &&(j<wskB)) 
				{
				PXX_mid+=2*FXP*X;
				PYY_mid+=2*FYP*Y;
				PZZ_mid+=2*FZP*Z; 
				PYZ_mid+=2*FYP*Z; 
				}
			PXY+=2*FXP*Y;
			PYX+=2*FYP*X;
			PXZ+=2*FXP*Z;  
			PZX+=2*FZP*X; 
			PZY+=2*FZP*Y;
			PYZ+=2*FYP*Z;	
		
			//------------------------------------------------------------------	 
			}   
			else if(R2>=RCUT*RCUT)
			{
			EUP=0.0; FXP=0.0; FYP=0.0; FZP=0.0;		 
			}   
			FX[i]=FX[i]+FXP; FY[i]=FY[i]+FYP;	FZ[i]=FZ[i]+FZP; U[i]+=EUP;
			FX[j]=FX[j]-FXP; FY[j]=FY[j]-FYP;	FZ[j]=FZ[j]-FZP; 
		//----------------------------------------------------------------------	 
			}
    		}
		if(fabs(RZ[i])>30)
		{
			FZ[i]=FZ[i]-100*(FZ[i]-30)*(FZ[i]-30);
		}
		
		
		}
		
		FABZ1=0;
		for(int ii=0; ii<wskA; ii++)
    	{
		FABZ1=FABZ1+FZ[ii];
		FABX=FABX+FX[ii];
		FABY=FABY+FY[ii];
		}
    	FABZ1=-2*FABZ1;
    	FABX=-2*FABX;
    	FABY=-2*FABY;
    	
    	FBCZ1=0;
		for(int ii=wskB; ii<NATOMS; ii++)
    	{
		FBCZ1=FBCZ1+FZ[ii];
		FBCX=FBCX+FX[ii];
		FBCY=FBCY+FY[ii];
		}
    	FBCZ1=2*FBCZ1;
    	FBCX=2*FBCX;
    	FBCY=2*FBCY;
    	
//------------------------------------------------------------------------------
sumPXX_mid=sumPXX_mid+PXX_mid;
sumPYY_mid=sumPYY_mid+PYY_mid;
sumPZZ_mid=sumPZZ_mid+PZZ_mid;
sumPYZ_mid=sumPYZ_mid+PYZ_mid;

meanPXX_mid=sumPXX_mid/DENOM;
meanPYY_mid=sumPYY_mid/DENOM;
meanPZZ_mid=sumPZZ_mid/DENOM;
meanPYZ_mid=sumPYZ_mid/DENOM;


//------------------------------------------------------------------------------
PVIR=pow(3*Volume,-1)*(0.5*(PXX+PYY+PZZ))+rho*T;

PVIR_mid=pow(3*Volum2,-1)*(0.5*(PXX_mid+PYY_mid+PZZ_mid))+(NATOMS*pow(Volum2,-1))*T_mid_2;
PVIR_mid2=pow(3*Volum4,-1)*(0.5*(PXX_mid+PYY_mid+PZZ_mid))+(NATOMS*pow(Volum4,-1))*T_mid_2;
PVIR_mid3=pow(3*Volum3,-1)*(0.5*(PXX_mid+PYY_mid+PZZ_mid))+(NATOMS*pow(Volum3,-1))*T_mid_2;

PVIR1=pow(3*Volum1,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(Volum1,-1))*T;
PVIR2=pow(3*Volum2,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(Volum2,-1))*T;
PVIR3=pow(3*Volum3,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(Volum3,-1))*T;
PVIR4=pow(3*Volum4,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(Volum4,-1))*T;

SPVIR1=pow(3*meanVolum1,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(meanVolum1,-1))*T;
SPVIR2=pow(3*meanVolum2,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(meanVolum2,-1))*T;
SPVIR3=pow(3*meanVolum3,-1)*(0.5*(PXX+PYY+PZZ))+(NATOMS*pow(meanVolum3,-1))*T;

sumPVIR1=sumPVIR1+PVIR1;
sumPVIR2=sumPVIR2+PVIR2;
sumPVIR3=sumPVIR3+PVIR3;
sumPVIR4=sumPVIR4+PVIR4;
sumSPVIR1=sumSPVIR1+SPVIR1;
sumSPVIR2=sumSPVIR2+SPVIR2;
sumSPVIR3=sumSPVIR3+SPVIR3;

sumPVIR_mid1=sumPVIR_mid1+PVIR_mid;
sumPVIR_mid2=sumPVIR_mid2+PVIR_mid2;
sumPVIR_mid3=sumPVIR_mid3+PVIR_mid3;

meanPVIR_mid1=sumPVIR_mid1*pow(DENOM,-1);
meanPVIR_mid2=sumPVIR_mid2*pow(DENOM,-1);
meanPVIR_mid3=sumPVIR_mid3*pow(DENOM,-1);

meanPVIR1=sumPVIR1*pow(DENOM,-1);
meanPVIR2=sumPVIR2*pow(DENOM,-1);
meanPVIR3=sumPVIR3*pow(DENOM,-1);
meanPVIR4=sumPVIR4*pow(DENOM,-1);
meanSPVIR1=sumSPVIR1*pow(DENOM,-1);
meanSPVIR2=sumSPVIR2*pow(DENOM,-1);
meanSPVIR3=sumSPVIR3*pow(DENOM,-1);

meanFTET1=sumFTET1*pow(DENOM,-1);
meanFTET2=sumFTET2*pow(DENOM,-1);

PVZ=pow(Volume,-1)*0.5*PZZ+EKZ;

sumPW=sumPW+PW;
meanPW=sumPW/DENOM;

sumPWX=sumPWX+PWX;
meanPWX=sumPWX/DENOM;

sumPWY=sumPWY+PWY;
meanPWY=sumPWY/DENOM;
//------------------------------------------------------------------------------
		for (int i=0;i<NATOMS;i++)
		{
		VXH1[i]=VXH2[i]; VYH1[i]=VYH2[i]; VZH1[i]=VZH2[i];
		RX_old[i]=RX[i]; RY_old[i]=RY[i]; RZ_old[i]=RZ[i];
		}
scale_all=1;
scale_up=1; scale_down=1;
		if(LOOP%freq_scale==0)
		{
			scale_down=T_down*pow(T0,-1);
			scale_up=T_up*pow(T0,-1);
			scale_all=T*pow(T0,-1);
		
		dEW1=(3.0/2.0)*wskA*T0*(scale_down-1);
		dEW2=(3.0/2.0)*wskA*T0*(scale_up-1);
		sumDEW1=sumDEW1+dEW1;
		sumDEW2=sumDEW2+dEW2;
		meanDEW1=sumDEW1*pow(DENOM,-1);
		meanDEW2=sumDEW2*pow(DENOM,-1);
		
			if(scale_vel_w==1)
			{
				for (int i=0; i<NATOMS; i++)
    			{  
    			VXH1[i]=VXH1[i]*pow(scale_all,-0.5); 
				VYH1[i]=VYH1[i]*pow(scale_all,-0.5); 
				VZH1[i]=VZH1[i]*pow(scale_all,-0.5);  
				}
			}
		
			else if(scale_vel_w==2)
			{
			//------------------------------------------------------------------	
				PPXA=0; PPYA=0; PPZA=0;
				for (int i=0; i<wskA; i++)
    			{  
    			PPXA=PPXA-VXH1[i];
    			PPYA=PPYA-VYH1[i];
    			PPZA=PPZA-VZH1[i];
    		
    			VXH1[i]=VXH1[i]*pow(scale_down,-0.5); 
				VYH1[i]=VYH1[i]*pow(scale_down,-0.5); 
				VZH1[i]=VZH1[i]*pow(scale_down,-0.5);  
				
				PPXA=PPXA+VXH1[i];
    			PPYA=PPYA+VYH1[i];
    			PPZA=PPZA+VZH1[i];
				}
				PPXA=PPXA/wskA;
				PPYA=PPYA/wskA;
				PPZA=PPZA/wskA;	
				
				for (int i=0; i<wskA; i++)
    			{      			
    			VXH1[i]=VXH1[i]-korekta_pedu_scian*PPXA; 
				VYH1[i]=VYH1[i]-korekta_pedu_scian*PPYA; 
				VZH1[i]=VZH1[i]-korekta_pedu_scian*PPZA;  
				}
			//------------------------------------------------------------------	
				
				PPXC=0; PPYC=0; PPZC=0;
				for (int i=wskB; i<NATOMS; i++)
    			{  
    			PPXC=PPXC-VXH1[i];
    			PPYC=PPYC-VYH1[i];
    			PPZC=PPZC-VZH1[i];
    			
    			VXH1[i]=VXH1[i]*pow(scale_up,-0.5); 
				VYH1[i]=VYH1[i]*pow(scale_up,-0.5); 
				VZH1[i]=VZH1[i]*pow(scale_up,-0.5);  
				
				PPXC=PPXC+VXH1[i];
    			PPYC=PPYC+VYH1[i];
    			PPZC=PPZC+VZH1[i];
				}		
				PPXC=PPXC/wskA;
				PPYC=PPYC/wskA;
				PPZC=PPZC/wskA;
				
				
				
				for (int i=wskB; i<NATOMS; i++)
    			{  
    			VXH1[i]=VXH1[i]-korekta_pedu_scian*PPXC; 
				VYH1[i]=VYH1[i]-korekta_pedu_scian*PPYC; 
				VZH1[i]=VZH1[i]-korekta_pedu_scian*PPZC;  
				}
			//------------------------------------------------------------------
				PPXB=0; PPYB=0; PPZB=0;
				for (int i=wskA; i<wskB; i++)
    			{
				PPXB=PPXB+VXH1[i];
    			PPYB=PPYB+VYH1[i];
    			PPZB=PPZB+VZH1[i];
				}
			
				PPXB=PPXB/(NATOMS-2*wskA);
				PPYB=PPYB/(NATOMS-2*wskA);
				PPZB=PPZB/(NATOMS-2*wskA);
			//------------------------------------------------------------------
			}
		
		} 

FEXT=P0*BX*BY/wskA;
FYext=FY0*BX*BY/wskA;


VA=0;VC=0;
VYA=0; VYC=0;
		for (int i=0;i<wskA;i++)
		{
			VA=VA+VZH1[i];
			VYA=VYA+VYH1[i];
		}
		VA=VA/wskA;
		VYA=VYA/wskA;
		
		for (int i=wskB;i<NATOMS;i++)
		{
			VC=VC+VZH1[i];
			VYC=VYC+VYH1[i];
		}
		VC=VC/wskA;
		VYC=VYC/wskA;


VW=0.5*(VYC-VYA);
sumVW=sumVW+VW;
meanVW=sumVW/DENOM;

XA=-(1/Qp)*(NEWP1-P0*BX*BY);

XC=-(1/Qp)*(NEWP2-P0*BX*BY);

ST=0;

distA=distA+VYA*DT;
distC=distC+VYC*DT;

shFA=0;
shFC=0;
//------------------------------------------------------------------------------
//WARSTWA A
//------------------------------------------------------------------------------
TERM=1;	
SHEAR=1;

		for (int i=0;i<wskA;i++)
		{
		
		VXH2[i]=VXH1[i]+DT*(FX[i]);
		VYH2[i]=VYH1[i]+DT*(FY[i]+SH1*gam2*(0.5*(VYC-VYA)-VSH));
		VZH2[i]=VZH1[i]+DT*(FZ[i]+BARO4*(FEXT-BARO5*beta*0.5*(VA-VC)));
		shFA=shFA+SH1*gam2*(0.5*(VYC-VYA)-VSH);
		
		ST=ST+SH1*0.5*P0*BX*BX*((VYC-VYA)-2*VSH)/(2*VSH);
	
		RY_eq[i]=RY_eq[i]+SHEAR*VSH*DT;
		RY_eq[i]=RY_eq[i]-YPERIOD*(DBOX*round(RY_eq[i]*pow(DBOX,-1)));
		
		RX[i]=RX_old[i]+DT*VXH2[i];
		RY[i]=RY_old[i]+DT*VYH2[i]+SH2*SHEAR*VSH*DT+SH3*(Qp*((VYC-VYA)-2*VSH))*0.5*DT*DT;
		RZ[i]=RZ_old[i]+DT*VZH2[i];
        
        DELTARX=0; DELTARY=0; DELTARZ=0;
        
        DELTARX=RX[i]-RX_old[i];
		DELTARY=RY[i]-RY_old[i];
		DELTARZ=RZ[i]-RZ_old[i];
		
		RXR[i]=RXR[i]+DELTARX;
		RYR[i]=RYR[i]+DELTARY;
		RZR[i]=RZR[i]+DELTARZ;
        
    	RX[i]=RX[i]-XPERIOD*(DBOX*round(RX[i]*pow(DBOX,-1)));
    	RY[i]=RY[i]-YPERIOD*(DBOX*round(RY[i]*pow(DBOX,-1)));
    	RZ[i]=RZ[i]-ZPERIOD*(DBOZ*round(RZ[i]*pow(DBOZ,-1)));
    	
    	//historgam pozycji Z---------------------------------------------------
    	WHZ=0;
    	WHZ=int(RX[i]*bin+0.5*((RX[i]+crc)*pow(fabs(RX[i]+crc),-1)))+0.5*RANG;
    	histX[WHZ]=histX[WHZ]+1*GET;
    	
    	WHZ=0;
        WHZ=int(RY[i]*bin+0.5*((RY[i]+crc)*pow(fabs(RY[i]+crc),-1)))+0.5*RANG;
    	histY[WHZ]=histY[WHZ]+1*GET;
    	
		WHZ=0;
    	WHZ=int(RZ[i]*bin+0.5*((RZ[i]+crc)*pow(fabs(RZ[i]+crc),-1)))+0.5*RANG;
    	histZ[WHZ]=histZ[WHZ]+1*GET;
    	//----------------------------------------------------------------------
	
		VX[i]=0.5*(VXH1[i]+VXH2[i]);
		VY[i]=0.5*(VYH1[i]+VYH2[i]);
		VZ[i]=0.5*(VZH1[i]+VZH2[i]);
		
		WHZ=0;
    	WHZ=int((int(RZ[i]*binRZ+0.5*((RZ[i]+crc)*pow(fabs(RZ[i]+crc),-1)))+0.5*TGAP)*pow(resol,-1)+0.5);
		
		VX_pec[i]=VX[i]-MEAN_RTVX[WHZ];
		VY_pec[i]=VY[i]-MEAN_RTVY[WHZ];
		VZ_pec[i]=VZ[i]-MEAN_RTVZ[WHZ];		
		
		VPEC[i]=pow(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i],0.5);
		
		WHZ=0;
    	WHZ=int(VX_pec[i]*binHV+0.5*((VX_pec[i]+crc)*pow(fabs(VX_pec[i]+crc),-1)))+0.5*TGAP;
    	histVXA[WHZ]=histVXA[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VY_pec[i]*binHV+0.5*((VY_pec[i]+crc)*pow(fabs(VY_pec[i]+crc),-1)))+0.5*TGAP;
    	histVYA[WHZ]=histVYA[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VZ_pec[i]*binHV+0.5*((VZ_pec[i]+crc)*pow(fabs(VZ_pec[i]+crc),-1)))+0.5*TGAP;
    	histVZA[WHZ]=histVZA[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VPEC[i]*binHV+0.5*((VPEC[i]+crc)*pow(fabs(VPEC[i]+crc),-1)))+0.5*TGAP;
    	histVA[WHZ]=histVA[WHZ]+1*GET;
		}
//------------------------------------------------------------------------------
//WARSTWA B
//------------------------------------------------------------------------------
TERM=1;	
SHEAR=0;
	
		for (int i=wskA;i<wskB;i++)
		{
		VXH2[i]=VXH1[i]+DT*(FX[i]);
		VYH2[i]=VYH1[i]+DT*(FY[i]);
		VZH2[i]=VZH1[i]+DT*(FZ[i]);
		
		RY_eq[i]=RY_eq[i]+SHEAR*VSH*DT;
		RY_eq[i]=RY_eq[i]-YPERIOD*(DBOX*round(RY_eq[i]*pow(DBOX,-1)));
		
		RX[i]=RX_old[i]+DT*VXH2[i];
		RY[i]=RY_old[i]+DT*VYH2[i];
		RZ[i]=RZ_old[i]+DT*VZH2[i];
        
        DELTARX=0; DELTARY=0; DELTARZ=0;
        
        DELTARX=RX[i]-RX_old[i];
		DELTARY=RY[i]-RY_old[i];
		DELTARZ=RZ[i]-RZ_old[i];
		
		RXR[i]=RXR[i]+DELTARX;
		RYR[i]=RYR[i]+DELTARY;
		RZR[i]=RZR[i]+DELTARZ;
        
    	RX[i]=RX[i]-XPERIOD*(DBOX*round(RX[i]*pow(DBOX,-1)));
    	RY[i]=RY[i]-YPERIOD*(DBOX*round(RY[i]*pow(DBOX,-1)));
    	RZ[i]=RZ[i]-ZPERIOD*(DBOZ*round(RZ[i]*pow(DBOZ,-1)));
    	
    	//historgam pozycji Z---------------------------------------------------
    	WHZ=0;
    	WHZ=int(RX[i]*bin+0.5*((RX[i]+crc)*pow(fabs(RX[i]+crc),-1)))+0.5*RANG;
    	histX[WHZ]=histX[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(RY[i]*bin+0.5*((RY[i]+crc)*pow(fabs(RY[i]+crc),-1)))+0.5*RANG;
    	histY[WHZ]=histY[WHZ]+1*GET;
    	
		WHZ=0;
    	WHZ=int(RZ[i]*bin+0.5*((RZ[i]+crc)*pow(fabs(RZ[i]+crc),-1)))+0.5*RANG;
    	histZ[WHZ]=histZ[WHZ]+1*GET;
    	histZ2[WHZ]=histZ2[WHZ]+1*GET;
    	//----------------------------------------------------------------------
	
		VX[i]=0.5*(VXH1[i]+VXH2[i]);
		VY[i]=0.5*(VYH1[i]+VYH2[i]);
		VZ[i]=0.5*(VZH1[i]+VZH2[i]);

		
		WHZ=0;
    	WHZ=int((int(RZ[i]*binRZ+0.5*((RZ[i]+crc)*pow(fabs(RZ[i]+crc),-1)))+0.5*TGAP)*pow(resol,-1)+0.5);
		
		VX_pec[i]=VX[i]-MEAN_RTVX[WHZ];
		VY_pec[i]=VY[i]-MEAN_RTVY[WHZ];
		VZ_pec[i]=VZ[i]-MEAN_RTVZ[WHZ];
		
		VPEC[i]=pow(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i],0.5);
		
		WHZ=0;
    	WHZ=int(VX_pec[i]*binHV+0.5*((VX_pec[i]+crc)*pow(fabs(VX_pec[i]+crc),-1)))+0.5*TGAP;
    	histVXB[WHZ]=histVXB[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VY_pec[i]*binHV+0.5*((VY_pec[i]+crc)*pow(fabs(VY_pec[i]+crc),-1)))+0.5*TGAP;
    	histVYB[WHZ]=histVYB[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VZ_pec[i]*binHV+0.5*((VZ_pec[i]+crc)*pow(fabs(VZ_pec[i]+crc),-1)))+0.5*TGAP;
    	histVZB[WHZ]=histVZB[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VPEC[i]*binHV+0.5*((VPEC[i]+crc)*pow(fabs(VPEC[i]+crc),-1)))+0.5*TGAP;
    	histVB[WHZ]=histVB[WHZ]+1*GET;
		}

//------------------------------------------------------------------------------
//WARSTWA C
//------------------------------------------------------------------------------		
TERM=1;	
SHEAR=-1;	

		for (int i=wskB;i<NATOMS;i++)
		{
		VXH2[i]=VXH1[i]+DT*(FX[i]);
		VYH2[i]=VYH1[i]+DT*(FY[i]-SH1*gam2*(0.5*(VYC-VYA)-VSH));
		VZH2[i]=VZH1[i]+DT*(FZ[i]-BARO4*(FEXT-BARO5*beta*0.5*(VA-VC)));
		
		shFC=shFC-SH1*gam2*(0.5*(VYC-VYA)-VSH);
		
		ST=ST-SH1*0.5*P0*BX*BX*((VYC-VYA)-2*VSH)/(2*VSH);		
		
		RY_eq[i]=RY_eq[i]+SHEAR*VSH*DT;
		RY_eq[i]=RY_eq[i]-YPERIOD*(DBOX*round(RY_eq[i]*pow(DBOX,-1)));
		
		RX[i]=RX_old[i]+DT*VXH2[i];
		RY[i]=RY_old[i]+DT*VYH2[i]+SH2*SHEAR*VSH*DT-SH3*(Qp*((VYC-VYA)-2*VSH))*0.5*DT*DT;
		RZ[i]=RZ_old[i]+DT*VZH2[i];
        
        DELTARX=0; DELTARY=0; DELTARZ=0;
        
        DELTARX=RX[i]-RX_old[i];
		DELTARY=RY[i]-RY_old[i];
		DELTARZ=RZ[i]-RZ_old[i];
		
		RXR[i]=RXR[i]+DELTARX;
		RYR[i]=RYR[i]+DELTARY;
		RZR[i]=RZR[i]+DELTARZ;
        
    	RX[i]=RX[i]-XPERIOD*(DBOX*round(RX[i]*pow(DBOX,-1)));
    	RY[i]=RY[i]-YPERIOD*(DBOX*round(RY[i]*pow(DBOX,-1)));
    	RZ[i]=RZ[i]-ZPERIOD*(DBOZ*round(RZ[i]*pow(DBOZ,-1)));
    	
    	//historgam pozycji Z---------------------------------------------------
    	WHZ=0;
    	WHZ=int(RX[i]*bin+0.5*((RX[i]+crc)*pow(fabs(RX[i]+crc),-1)))+0.5*RANG;
    	histX[WHZ]=histX[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(RY[i]*bin+0.5*((RY[i]+crc)*pow(fabs(RY[i]+crc),-1)))+0.5*RANG;
    	histY[WHZ]=histY[WHZ]+1*GET;
    	
		WHZ=0;
    	WHZ=int(RZ[i]*bin+0.5*((RZ[i]+crc)*pow(fabs(RZ[i]+crc),-1)))+0.5*RANG;
		histZ[WHZ]=histZ[WHZ]+1*GET;
    	//----------------------------------------------------------------------
	
		VX[i]=0.5*(VXH1[i]+VXH2[i]);
		VY[i]=0.5*(VYH1[i]+VYH2[i]);
		VZ[i]=0.5*(VZH1[i]+VZH2[i]);
		
		WHZ=0;
    	WHZ=int((int(RZ[i]*binRZ+0.5*((RZ[i]+crc)*pow(fabs(RZ[i]+crc),-1)))+0.5*TGAP)*pow(resol,-1)+0.5);
		
		VX_pec[i]=VX[i]-MEAN_RTVX[WHZ];
		VY_pec[i]=VY[i]-MEAN_RTVY[WHZ];
		VZ_pec[i]=VZ[i]-MEAN_RTVZ[WHZ];		
		
		VPEC[i]=pow(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i],0.5);
		
		WHZ=0;
    	WHZ=int(VX_pec[i]*binHV+0.5*((VX_pec[i]+crc)*pow(fabs(VX_pec[i]+crc),-1)))+0.5*TGAP;
    	histVXC[WHZ]=histVXC[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VY_pec[i]*binHV+0.5*((VY_pec[i]+crc)*pow(fabs(VY_pec[i]+crc),-1)))+0.5*TGAP;
    	histVYC[WHZ]=histVYC[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VZ_pec[i]*binHV+0.5*((VZ_pec[i]+crc)*pow(fabs(VZ_pec[i]+crc),-1)))+0.5*TGAP;
    	histVZC[WHZ]=histVZC[WHZ]+1*GET;
    	
    	WHZ=0;
    	WHZ=int(VPEC[i]*binHV+0.5*((VPEC[i]+crc)*pow(fabs(VPEC[i]+crc),-1)))+0.5*TGAP;
    	histVC[WHZ]=histVC[WHZ]+1*GET;
		
		}
//------------------------------------------------------------------------------
SXP=0; SYP=0; SZP=0;
		for (int i=0; i<NATOMS; i++) 
        {    		
    		SXP=SXP+RX[i];
    		SYP=SYP+RY[i];
    		SZP=SZP+RZ[i];
        }
SXP=SXP*pow(NATOMS,-1); SYP=SYP*pow(NATOMS,-1); SZP=SZP*pow(NATOMS,-1);     

if(zerowanieSM2==1)
{
for (int i=0; i<NATOMS; i++) 
        {
            RX[i]=RX[i]-SXP;
            RY[i]=RY[i]-SYP;
            RZ[i]=RZ[i]-SZP;
            
            RX_eq[i]=RX[i];
    		RY_eq[i]=RY[i];
  		    RZ_eq[i]=RZ[i];
   		
       }
}
SXP=0; SYP=0; SZP=0;
	
//TEMPERATURA W WYSOKOSCI-------------------------------------------------------	
for(int j=0;j<TGAP;j++)
{
	
RTZ[j]=0;
COUNT[j]=0;

RTZ_2[j]=0;
COUNT_2[j]=0;

RTZ_X[j]=0;
COUNT_X[j]=0;
RTZ_Y[j]=0;
COUNT_Y[j]=0;
RTZ_Z[j]=0;
COUNT_Z[j]=0;

RVX[j]=0;
C_VX[j]=0;

RVY[j]=0;
C_VY[j]=0;

RVY_pec[j]=0;
C_VY_pec[j]=0;

RVZ[j]=0;
C_VZ[j]=0;

	for(int i=0;i<NATOMS;i++)
	{
	double indeks;
	indeks=(-0.5*TGAP+j)*resol;	
	if(binRZ*RZ[i]>indeks && binRZ*RZ[i]<indeks+resol)
	{	
	RTZ[j]=RTZ[j]+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);
	COUNT[j]=COUNT[j]+1;
	
	RTZ_2[j]=RTZ_2[j]+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);
	COUNT_2[j]=COUNT_2[j]+1;
	
	RTZ_X[j]=RTZ_X[j]+0.5*(VX[i]*VX[i]);
	COUNT_X[j]=COUNT_X[j]+1;
	
	RTZ_Y[j]=RTZ_Y[j]+0.5*(VY[i]*VY[i]);
	COUNT_Y[j]=COUNT_Y[j]+1;
	
	RTZ_Z[j]=RTZ_Z[j]+0.5*(VZ[i]*VZ[i]);
	COUNT_Z[j]=COUNT_Z[j]+1;
	
	RVX[j]=RVX[j]+VX[i];
	C_VX[j]=C_VX[j]+1;
	
	RVY[j]=RVY[j]+VY[i];
	C_VY[j]=C_VY[j]+1;
	
	RVY_pec[j]=RVY_pec[j]+VY_pec[i];
	C_VY_pec[j]=C_VY_pec[j]+1;
	
	RVZ[j]=RVZ[j]+VZ[i];
	C_VZ[j]=C_VZ[j]+1;
	}

	}
//------------------------------------------------------------------------------	
if(COUNT[j]>0)
{	
	RTZ[j]=RTZ[j]*2*pow(3*COUNT[j],-1);
	preTabDenomRTZ[j]=1;
}
else
{
	RTZ[j]=0;
	preTabDenomRTZ[j]=0;
}
SUM_RTZ[j]=SUM_RTZ[j]+RTZ[j];
MEAN_RTZ[j]=SUM_RTZ[j]*pow(DENOM,-1);	
tabDenomRTZ[j]=tabDenomRTZ[j]+preTabDenomRTZ[j];
//------------------------------------------------------------------------------	
if(COUNT_2[j]>0)
{	
	RTZ_2[j]=RTZ_2[j]*2*pow(3*COUNT_2[j],-1);
	preTabDenomRTZ_2[j]=1;
}
else
{
	RTZ_2[j]=0;
	preTabDenomRTZ_2[j]=0;
}
SUM_RTZ_2[j]=SUM_RTZ_2[j]+RTZ_2[j];
MEAN_RTZ_2[j]=SUM_RTZ_2[j]*pow(DENOM,-1);	
tabDenomRTZ_2[j]=tabDenomRTZ_2[j]+preTabDenomRTZ_2[j]; 
//------------------------------------------------------------------------------
if(COUNT_X[j]>0)
{	
	RTZ_X[j]=RTZ_X[j]*2*pow(3*COUNT_X[j],-1);
	preTabDenomRTZ_X[j]=1;
}
else
{
	RTZ_X[j]=0;
	preTabDenomRTZ_X[j]=0;
}
SUM_RTZ_X[j]=SUM_RTZ_X[j]+RTZ_X[j];
MEAN_RTZ_X[j]=SUM_RTZ_X[j]*pow(DENOM,-1);	
tabDenomRTZ_X[j]=tabDenomRTZ_X[j]+preTabDenomRTZ_X[j];
//------------------------------------------------------------------------------
if(COUNT_Y[j]>0)
{	
	RTZ_Y[j]=RTZ_Y[j]*2*pow(3*COUNT_Y[j],-1);
	preTabDenomRTZ_Y[j]=1;
}
else
{
	RTZ_Y[j]=0;
	preTabDenomRTZ_Y[j]=0;
}
SUM_RTZ_Y[j]=SUM_RTZ_Y[j]+RTZ_Y[j];
MEAN_RTZ_Y[j]=SUM_RTZ_Y[j]*pow(DENOM,-1);
tabDenomRTZ_Y[j]=tabDenomRTZ_Y[j]+preTabDenomRTZ_Y[j]; 
//------------------------------------------------------------------------------
if(COUNT_Z[j]>0)
{	
	RTZ_Z[j]=RTZ_Z[j]*2*pow(3*COUNT_Z[j],-1);
	preTabDenomRTZ_Z[j]=1;
}
else
{
	RTZ_Z[j]=0;
	preTabDenomRTZ_Z[j]=0;
}
SUM_RTZ_Z[j]=SUM_RTZ_Z[j]+RTZ_Z[j];
MEAN_RTZ_Z[j]=SUM_RTZ_Z[j]*pow(DENOM,-1);
tabDenomRTZ_Z[j]=tabDenomRTZ_Z[j]+preTabDenomRTZ_Z[j];
//------------------------------------------------------------------------------
if(C_VX[j]>0)
{	
	RVX[j]=RVX[j]*pow(C_VX[j],-1);
	preTabDenomRVX[j]=1;
}
else
{
	RVX[j]=0;
	preTabDenomRVX[j]=1;
}
SUM_RTVX[j]=SUM_RTVX[j]+RVX[j], 
MEAN_RTVX[j]=SUM_RTVX[j]*pow(DENOM,-1);
tabDenomRVX[j]=tabDenomRVX[j]+preTabDenomRVX[j];
//------------------------------------------------------------------------------
if(C_VY[j]>0)
{	
	RVY[j]=RVY[j]*pow(C_VY[j],-1);
	preTabDenomRVY[j]=1;
}
else
{
	RVY[j]=0;
	preTabDenomRVY[j]=0;
}
SUM_RTVY[j]=SUM_RTVY[j]+RVY[j], 
MEAN_RTVY[j]=SUM_RTVY[j]*pow(DENOM,-1);
tabDenomRVY[j]=tabDenomRVY[j]+preTabDenomRVY[j];
//------------------------------------------------------------------------------
if(C_VY_pec[j]>0)
{	
	RVY_pec[j]=RVY_pec[j]*pow(C_VY_pec[j],-1);
	preTabDenomRVY_pec[j]=1;
}
else
{
	RVY_pec[j]=0;
	preTabDenomRVY_pec[j]=0;
}
SUM_RTVY_pec[j]=SUM_RTVY_pec[j]+RVY_pec[j], 
MEAN_RTVY_pec[j]=SUM_RTVY_pec[j]*pow(n+1,-1);
tabDenomRVY_pec[j]=tabDenomRVY_pec[j]+preTabDenomRVY_pec[j]; 
//------------------------------------------------------------------------------
if(C_VZ[j]>0)
{	
	RVZ[j]=RVZ[j]*pow(C_VZ[j],-1);
	preTabDenomRVZ[j]=1;
}
else
{
	RVZ[j]=0;
	preTabDenomRVZ[j]=0;
}
SUM_RTVZ[j]=SUM_RTVZ[j]+RVZ[j], 
MEAN_RTVZ[j]=SUM_RTVZ[j]*pow(DENOM,-1);
tabDenomRVZ[j]=tabDenomRVZ[j]+preTabDenomRVZ[j]; 
//------------------------------------------------------------------------------
}	
	
//BAROSTAT----------------------------------------------------------------------

F=1*pow(1+DT*pow(tauH,-1),-1);
E=F*(1-DT*pow(tauH,-1));
MA=200;

A=BX*BX;

NEWP1=0.5*(FABZ1)*pow(DBOX*DBOX,-1);
NEWP2=0.5*(FBCZ1)*pow(DBOX*DBOX,-1);

NEWPX1=0.5*(FABX)*pow(DBOX*DBOX,-1);
NEWPX2=0.5*(FBCX)*pow(DBOX*DBOX,-1);

NEWPY1=0.5*(FABY)*pow(DBOX*DBOX,-1);
NEWPY2=0.5*(FBCY)*pow(DBOX*DBOX,-1);

PW=0.5*(NEWP1+NEWP2);

if(cisnienieNaturalne==1 && n==0)
{
	P0=PW;
	P0=round(100*P0)/100;
	cout<<"P0 przestawione na wartosc: "<<P0<<endl;
}

PWX=0.5*(NEWPX1+NEWPX2);
PWY=0.5*(NEWPY1+NEWPY2);

	Pt=(FUP2-FDOWN2)*pow((DBOX*DBOX),-1);
	Pt=0.5*Pt;
	PTX=(FXRIGHT-FXLEFT)*pow((DBOX*DBOX),-1);
	PTY=(FYRIGHT-FYLEFT)*pow((DBOX*DBOX),-1);
	PTZ=Pt;
	
	
//KINETYCZNY WKLAD MOPu---------------------------------------------------------	
sumk=0;
if(sumk==0)
{
liczsumk=liczsumk+1;
}
	
//Vol av------------------------------------------------------------------------	
VAPstart=0;
VAPend=PGAP;
VAP=0;
VAPcount=0;
VAP2=0;
VAP2count=0;
VAP3=0;
VAP3count=0;
VAP4=0;
VAP4count=0;
VAP5=0;
VAP5count=0;
VAP10=0;
VAP10count=0;
VAP20=0;
VAP20count=0;
for(int iii=VAPstart; iii<VAPend; iii++)
{
	VAP=VAP+PMOP[iii];
	VAPcount=VAPcount+1; 	
	
	if((iii+1)%2==0)
	{
	VAP2=VAP2+PMOP[iii];
	VAP2count=VAP2count+1; 	
	}
	if((iii+1)%3==0)
	{
	VAP3=VAP3+PMOP[iii];
	VAP3count=VAP3count+1; 	
	}
	if((iii+1)%4==0)
	{
	VAP4=VAP4+PMOP[iii];
	VAP4count=VAP4count+1; 	
	}
	if((iii+1)%5==0)
	{
	VAP5=VAP5+PMOP[iii];
	VAP5count=VAP5count+1; 	
	}
	if((iii+1)%10==0)
	{
	VAP10=VAP10+PMOP[iii];
	VAP10count=VAP10count+1; 	
	}
	if(iii==0.5*PGAP)
	{
	VAP20=VAP20+PMOP[iii];
	VAP20count=VAP20count+1; 	
	}
	
}	

VAP=VAP*pow(VAPcount,-1);

VAP2=VAP2*pow(VAP2count,-1);
VAP3=VAP3*pow(VAP3count,-1);
VAP4=VAP4*pow(VAP4count,-1);
VAP5=VAP5*pow(VAP5count,-1);
VAP10=VAP10*pow(VAP10count,-1);
VAP20=VAP20*pow(VAP20count,-1);

	
//------------------------------------------------------------------------------	
	sumP=sumP+Pt;
	
	sumNEWP1=sumNEWP1+NEWP1;
	sumNEWP2=sumNEWP2+NEWP2;
	
	meanNEWP1=sumNEWP1*pow(DENOM,-1);
	meanNEWP2=sumNEWP2*pow(DENOM,-1);

//------------------------------------------------------------------------------	
	vzz=0;
	vzz_dec=BAR*vzz;

	V_down=0; V_up=0;
	for(int i=0;i<wskA;i++)
	{
	RZ[i]=RZ[i]-0.5*vzz_dec*DT;
	RZ_eq[i]=RZ_eq[i]-0.5*vzz_dec*DT*CONSTR;
	V_down=V_down+VZ[i];
	}

	for(int i=wskB;i<NATOMS;i++)
	{
	RZ[i]=RZ[i]+0.5*vzz_dec*DT;
	RZ_eq[i]=RZ_eq[i]+0.5*vzz_dec*DT*CONSTR;
	V_up=V_up+VZ[i];
	}
//------------------------------------------------------------------------------
	if(CONSTR==2)
	{
		VWD=FABZ1*pow(wskA,-1);
		for(int i=0;i<wskA;i++)
		{
		RZ_eq[i]=RZ_eq[i]+0.5*VWD*DT;
		}
		VWU=FBCZ1*pow(wskA,-1);
		for(int i=wskB;i<NATOMS;i++)
		{
		RZ_eq[i]=RZ_eq[i]-0.5*VWD*DT;
		}
	}
	
//aktualizacja listy Verleta----------------------------------------------------	

for(int i=0;i<NATOMS;i++)
	{	
		DRX[i]=DRX[i]+(RX[i]-RX_old[i]);
		DRY[i]=DRY[i]+(RY[i]-RY_old[i]);
		DRZ[i]=DRZ[i]+(RZ[i]-RZ_old[i]);
	
		DRX[i]=DRX[i]-XPERIOD*(DBOX*round(DRX[i]*pow(DBOX,-1))); 
		DRY[i]=DRY[i]-YPERIOD*(DBOX*round(DRY[i]*pow(DBOX,-1))); 
		DRZ[i]=DRZ[i]-ZPERIOD*(DBOZ*round(DRZ[i]*pow(DBOZ,-1)));
		
		if(DRX[i]*DRX[i]+DRY[i]*DRY[i]+DRZ[i]*DRZ[i]>(0.5*SKIN)*(0.5*SKIN))
		{
			counter=1;
			cout<<"Aktualizacja listy Verleta..."<<endl;
			cout<<"n: "<<n<<"  i: "<<i<<"  counter: "<<counter<<"  DR: "<<pow(DRX[i]*DRX[i]+DRY[i]*DRY[i]+DRZ[i]*DRZ[i],0.5)<<endl;
		}
	}

//------------------------------------------------------------------------------
ELJ=0;
EPH=0;
for(int i=0;i<NATOMS;i++)
{
ELJ=ELJ+ULJ[i];
EPH=EPH+UPH[i];	
}	
	
	EPC_walls=0;
	for(int i=0;i<wskA;i++) {EPC_walls=EPC_walls+U[i];}
	for(int i=wskB;i<NATOMS;i++) {EPC_walls=EPC_walls+U[i];}
	
	EPC_w_down=0;
	for(int i=0;i<wskA;i++) {EPC_w_down=EPC_w_down+U[i];}

	EPC_w_up=0;
	for(int i=wskB;i<NATOMS;i++) {EPC_w_up=EPC_w_up+U[i];}

	EPC_inter=0;
	for(int i=wskA;i<wskB;i++) {EPC_inter=EPC_inter+U[i];}
	
	EPC=0;
	for(int i=0;i<NATOMS;i++) {EPC=EPC+U[i];}

	EKC=0;
	g_in=0;
	for(int i=0;i<wskA;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
	g_in=3*(wskA);
	T_down=2*pow(g_in,-1)*EKC;

	EKC=0;
	g_in=0;
	for(int i=wskB;i<NATOMS;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
	g_in=3*(wskA);
	T_up=2*pow(g_in,-1)*EKC;

	EKC=0;
	g_in=0;
	for(int i=wskA;i<wskB;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	g_in=3*(NATOMS-2*wskA);
	T_mid=2*pow(g_in,-1)*EKC;
	
	EKC=0;
	g_in=0;
	for(int i=wskA;i<wskB;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
	g_in=3*(NATOMS-2*wskA);
	T_mid_2=2*pow(g_in,-1)*EKC;
	//-
	EKC=0;
	g_in=NATOMS;
	for(int i=0;i<wskA;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
	for(int i=wskA;i<wskB;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
	for(int i=wskB;i<NATOMS;i++) {EKC=EKC+0.5*(VX_pec[i]*VX_pec[i]+VY_pec[i]*VY_pec[i]+VZ_pec[i]*VZ_pec[i]);}
	T_pec=2*pow(g_in,-1)*EKC;
	//-
	EKC=0;
	g_in=0;
	for(int i=0;i<wskA;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	for(int i=wskB;i<NATOMS;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	g_in=3*(2*wskA);
	T_walls=2*pow(g_in,-1)*EKC;
//------------------------------------------------------------------------------	
	
	ENKX=0;
	g_in=3*NATOMS;
	for(int i=0;i<NATOMS;i++) {ENKX=ENKX+0.5*(VX[i]*VX[i]);}
	TEMPX=2*pow(g_in,-1)*ENKX;
	
	ENKY=0;
	g_in=3*NATOMS;
	for(int i=0;i<NATOMS;i++) {ENKY=ENKY+0.5*(VY[i]*VY[i]);}
	TEMPY=2*pow(g_in,-1)*ENKY;
	
	ENKZ=0;
	g_in=3*NATOMS;
	for(int i=0;i<NATOMS;i++) {ENKZ=ENKZ+0.5*(VZ[i]*VZ[i]);}
	TEMPZ=2*pow(g_in,-1)*ENKZ;
	
//------------------------------------------------------------------------------
	EKC=0;
	for(int i=0;i<NATOMS;i++) {EKC=EKC+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	T=2*pow(g,-1)*EKC;
	
	EKZ=0;
	for(int i=0;i<NATOMS;i++) {EKZ=EKZ+0.5*(VZ[i]*VZ[i]);}
	EKZ=EKZ*pow(NATOMS,-1);
	EKZ=EKZ*2*NATOMS*pow(Volume,-1);
	
	EKC_walls=0;
	for(int i=0;i<wskA;i++) {EKC_walls=EKC_walls+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	for(int i=wskB;i<NATOMS;i++) {EKC_walls=EKC_walls+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	
	EKC_w_down=0;
	for(int i=0;i<wskA;i++) {EKC_w_down=EKC_w_down+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	
	EKC_w_up=0;
	for(int i=wskB;i<NATOMS;i++) {EKC_w_up=EKC_w_up+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}

	EKC_inter=0;
	for(int i=wskA;i<wskB;i++) {EKC_inter=EKC_inter+0.5*(VX[i]*VX[i]+VY[i]*VY[i]+VZ[i]*VZ[i]);}
	
	calka=(EPC+EKC)*pow(NATOMS,-1);
	DH=pow(calkaBEGIN-calka,2);
	sumDH=sumDH+DH;
	meanDH=pow(sumDH*pow(DENOM,-1),0.5);


//mean values-------------------------------------------------------------------

SUM_EP=SUM_EP+EPC;
SUM_T=SUM_T+T;
SUM_T_up=SUM_T_up+T_up;
SUM_T_down=SUM_T_down+T_down;
SUM_T_mid=SUM_T_mid+T_mid;
SUM_T_mid_2=SUM_T_mid_2+T_mid_2;
SUM_T_walls=SUM_T_walls+T_walls;
SUM_T_pec=SUM_T_pec+T_pec;

sumTX=sumTX+TEMPX;
sumTY=sumTY+TEMPY;
sumTZ=sumTZ+TEMPZ;

SUM_P=SUM_P+PVIR;	

sumEPC_inter=sumEPC_inter+EPC_inter;
sumEPC_w_down=sumEPC_w_down+EPC_w_down;
sumEPC_w_up=sumEPC_w_up+EPC_w_up;

sumVAP=sumVAP+VAP;
sumVAP2=sumVAP2+VAP2;
sumVAP3=sumVAP3+VAP3;
sumVAP4=sumVAP4+VAP4;
sumVAP5=sumVAP5+VAP5;
sumVAP10=sumVAP10+VAP10;
sumVAP20=sumVAP20+VAP20;


meanEP=SUM_EP*pow(DENOM,-1);
meanT=SUM_T*pow(DENOM,-1);
meanT_up=SUM_T_up*pow(DENOM,-1);
meanT_down=SUM_T_down*pow(DENOM,-1);
meanT_mid=SUM_T_mid_2*pow(DENOM,-1);
meanT_mid_2=SUM_T_mid_2*pow(DENOM,-1);
meanT_walls=SUM_T_walls*pow(DENOM,-1);
meanP=SUM_P*pow(DENOM,-1);
meanT_pec=SUM_T_pec*pow(DENOM,-1);

meanTX=sumTX*pow(DENOM,-1);
meanTY=sumTY*pow(DENOM,-1);
meanTZ=sumTZ*pow(DENOM,-1);

meanEPC_inter=sumEPC_inter*pow(DENOM,-1); 
meanEPC_w_down=sumEPC_w_down*pow(DENOM,-1);
meanEPC_w_up=sumEPC_w_up*pow(DENOM,-1);

meanDR=meanDR*pow(DENOM,-1)*pow(NATOMS,-1);

meanVAP=sumVAP*pow(DENOM,-1);
meanVAP2=sumVAP2*pow(DENOM,-1);
meanVAP3=sumVAP3*pow(DENOM,-1);
meanVAP4=sumVAP4*pow(DENOM,-1);
meanVAP5=sumVAP5*pow(DENOM,-1);
meanVAP10=sumVAP10*pow(DENOM,-1);
meanVAP20=sumVAP20*pow(DENOM,-1);
//------------------------------------------------------------------------------

sumPXX=sumPXX+PXX;
sumPYY=sumPYY+PYY;
sumPZZ=sumPZZ+PZZ;
sumPZY=sumPZY+PZY;
sumPYZ=sumPYZ+PYZ;
sumPXY=sumPXY+PXY;
sumPYX=sumPYX+PYX;
sumPXZ=sumPXZ+PXZ;
sumPZX=sumPZX+PZX;

meanPXX=sumPXX/DENOM;
meanPYY=sumPYY/DENOM;
meanPZZ=sumPZZ/DENOM;
meanPZY=sumPZY/DENOM;
meanPYZ=sumPYZ/DENOM;
meanPXY=sumPXY/DENOM;
meanPYX=sumPYX/DENOM;
meanPXZ=sumPXZ/DENOM;
meanPZX=sumPZX/DENOM;

sumPTX=sumPTX+PTX;
sumPTY=sumPTY+PTY;

sumPTZ=sumPTZ+PTZ;
sumPVZ=sumPVZ+PVZ;
sumDIM=sumDIM+DIMZ;

meanPt=sumP*pow(DENOM,-1);

meanPZ=sumPTZ*pow(DENOM,-1);

meanL=sumL/DENOM;

SUM_PXX=SUM_PXX+PXX;
SUM_PYY=SUM_PYY+PYY;
SUM_PZZ=SUM_PZZ+PZZ;
SUM_PXY=SUM_PXY+PXY;
SUM_PYX=SUM_PYX+PYX;
SUM_PZX=SUM_PZX+PZX;
SUM_PXZ=SUM_PXZ+PXZ;
SUM_PYZ=SUM_PYZ+PYZ;
SUM_PZY=SUM_PZY+PZY;
SUM_FDOWN2=SUM_FDOWN2+FDOWN2;
SUM_FUP2=SUM_FUP2+FUP2;
SUM_FBXD=SUM_FBXD+FBXD;
SUM_FBXU=SUM_FBXU+FBXU;
SUM_FBYD=SUM_FBYD+FBYD;
SUM_FBYU=SUM_FBYU+FBYU;

PT1=FTET1*pow(A*wskA,-1);
PT2=FTET2*pow(A*wskA,-1); 
sumPT1=sumPT1+PT1;
sumPT2=sumPT2+PT2;

meanPT1=sumPT1*pow(DENOM,-1);
meanPT2=sumPT2*pow(DENOM,-1);

//constant of motions-----------------------------------------------------------
CPX=0; CPY=0; CPZ=0; CPXA=0; CPYA=0; CPZA=0; CPXB=0; CPYB=0; CPZB=0;
CPXC=0; CPYC=0; CPZC=0; SMX=0; SMY=0; SMZ=0; SMXA=0; SMYA=0; SMZA=0; 
SMXB=0; SMYB=0; SMZB=0; SMXC=0; SMYC=0; SMZC=0; CMPX=0; CMPY=0; CMPZ=0;

for (int ii=0;ii<NATOMS;ii++)
{
	if(ii<wskA)
	{
	CPXA=CPXA+VX[ii]; CPYA=CPYA+VY[ii]; CPZA=CPZA+VZ[ii];	
	SMXA=SMXA+RX[ii]; SMYA=SMYA+RY[ii]; SMZA=SMZA+RZ[ii];
	}
	else if(ii>=wskB)
	{
	CPXC=CPXC+VX[ii]; CPYC=CPYC+VY[ii]; CPZC=CPZC+VZ[ii];	
	SMXC=SMXC+RX[ii]; SMYC=SMYC+RY[ii]; SMZC=SMZC+RZ[ii];
	}
	else
	{
	CPXB=CPXB+VX[ii]; CPYB=CPYB+VY[ii]; CPZB=CPZB+VZ[ii];
	SMXB=SMXB+RX[ii]; SMYB=SMYB+RY[ii]; SMZB=SMZB+RZ[ii];	
	}
	
CPX=CPX+VX[ii]; CPY=CPY+VY[ii]; CPZ=CPZ+VZ[ii];
SMX=SMX+RX[ii]; SMY=SMY+RY[ii]; SMZ=SMZ+RZ[ii];

CMPX=CMPX+(VZ[ii]*RY[ii]-VY[ii]*RZ[ii]);
CMPY=CMPY+(VX[ii]*RZ[ii]-VZ[ii]*RX[ii]);
CMPZ=CMPZ+(VY[ii]*RX[ii]-VX[ii]*RY[ii]);
}

SMX=SMX*pow(NATOMS,-1);
SMY=SMY*pow(NATOMS,-1);
SMZ=SMZ*pow(NATOMS,-1);

SMXA=SMXA*pow(NATOMS,-1);
SMYA=SMYA*pow(NATOMS,-1);
SMZA=SMZA*pow(NATOMS,-1);

SMXB=SMXB*pow(NATOMS,-1);
SMYB=SMYB*pow(NATOMS,-1);
SMZB=SMZB*pow(NATOMS,-1);

SMXC=SMXC*pow(NATOMS,-1);
SMYC=SMYC*pow(NATOMS,-1);
SMZC=SMZC*pow(NATOMS,-1);

sumSMX=sumSMX+SMX;
sumSMY=sumSMY+SMY;
sumSMZ=sumSMZ+SMZ;
sumSMXA=sumSMXA+SMXA;
sumSMYA=sumSMYA+SMYA; 
sumSMZA=sumSMZA+SMZA; 
sumSMXB=sumSMXB+SMXB;
sumSMYB=sumSMYB+SMYB; 
sumSMZB=sumSMZB+SMZB; 
sumSMXC=sumSMXC+SMXC; 
sumSMYC=sumSMYC+SMYC; 
sumSMZC=sumSMZC+SMZC;

meanSMX=sumSMX*pow(DENOM,-1);
meanSMY=sumSMY*pow(DENOM,-1);
meanSMZ=sumSMZ*pow(DENOM,-1);
meanSMXA=sumSMXA*pow(DENOM,-1);
meanSMYA=sumSMYA*pow(DENOM,-1); 
meanSMZA=sumSMZA*pow(DENOM,-1); 
meanSMXB=sumSMXB*pow(DENOM,-1);
meanSMYB=sumSMYB*pow(DENOM,-1); 
meanSMZB=sumSMZB*pow(DENOM,-1); 
meanSMXC=sumSMXC*pow(DENOM,-1); 
meanSMYC=sumSMYC*pow(DENOM,-1); 
meanSMZC=sumSMZC*pow(DENOM,-1);
//------------------------------------------------------------------------------	
sumShFA=sumShFA+shFA;
sumShFC=sumShFC+shFC;

meanShFA=sumShFA/DENOM;
meanShFC=sumShFC/DENOM;	
//friction coef-----------------------------------------------------------------

MI1=-PWY/PW;
sumMI1=sumMI1+MI1;
meanMI1=sumMI1/DENOM;

MI2=-meanPZY*pow(meanPZZ,-1);
sumMI2=sumMI2+MI2;
meanMI2=sumMI2/DENOM;

MI3=-meanPYZ_mid*pow(meanPZZ_mid,-1);
sumMI3=sumMI3+MI3;
meanMI3=sumMI3/DENOM;

MI4=0.5*(shFC/(PW*BX*BX)-shFA/(PW*BX*BX));
sumMI4=sumMI4+MI4;
meanMI4=sumMI4/DENOM;

MI5=(meanDEW1+meanDEW2)/(BX*BX*PW*((VYC)-(VYA))*DT);
sumMI5=sumMI5+MI5;
meanMI5=sumMI5/DENOM;

SUMFRIC1=SUMFRIC1+MI1;
SUMFRIC2=SUMFRIC2+MI2; 
SUMFRIC3=SUMFRIC3+MI3;
SUMFRIC4=SUMFRIC4+MI4;

MI=-wskA*gam2*(0.5*(VYC-VYA)-VSH)/(PW*BX*BX);
sumMI=sumMI+MI;
meanMI=sumMI/DENOM;


//HISTOGRAMY TERMODYNAMICZNE
//------------------------------------------------------------------------------
WH=0;
WH=int(binT*T+0.5*((T+crc)*pow(fabs(T+crc),-1)));
histT[WH]=histT[WH]+GET*1;

WH=0; 
WH=int(binT*T_up+0.5*((T_up+crc)*pow(fabs(T_up+crc),-1)));
histT_up[WH]=histT_up[WH]+GET*1;

WH=0; 
WH=int(binT*T_down+0.5*((T_down+crc)*pow(fabs(T_down+crc),-1)));
histT_down[WH]=histT_down[WH]+GET*1; 

WH=0; 
WH=int(binT*T_mid+0.5*((T_mid+crc)*pow(fabs(T_mid+crc),-1)));
histT_mid[WH]=histT_mid[WH]+GET*1;

WH=0; 
WH=int(binT*T_mid_2+0.5*((T_mid_2+crc)*pow(fabs(T_mid_2+crc),-1)));
histT_mid_2[WH]=histT_mid_2[WH]+GET*1;

WH=0; 
WH=int(binL*(LENGTH-L0)+0.5*(((LENGTH-L0)+crc)*pow(fabs((LENGTH-L0)+crc),-1)));
histL[WH]=histL[WH]+GET*1;

WH=0; 
WH=int(binP*(PW-Pb)+0.5*(((PW-Pb)+crc)*pow(fabs((PW-Pb)+crc),-1)));
histP[WH]=histP[WH]+GET*1;

CNT=CNT+GET*1;

//ZAPIS DANYCH------------------------------------------------------------------
LOOP=n-1;
if(LOOP%JUMP2==0)
{
//------------------------------------------------------------------------------        
przebiegi.precision(16);
przebiegi<<
/*1*/n<<				"  "		<</*2*/tim<<					"  "		<</*3*/T<<"  "<<
/*4*/meanT<<			"  "		<</*5*/PW<<						"  "		<</*6*/meanPW<<"  "<<
/*7*/LENGTH<<			"  "		<</*8*/meanL<<					"  "		<</*9*/MI<<"  "<<
/*10*/meanMI<<			"  "		<</*11*/0.5*(VYC-VYA)<<			"  "		<</*12*/BX*BY<<"\n";
przebiegi.flush();

cout<<"Wykonano: "<<n<<" krokow czasowych z: "<<RUN<< " wszystkich krokow..."<<endl;

//------------------------------------------------------------------------------        
		}	
		

//ZAPIS POZYCJI I PROFILI-------------------------------------------------------

        if(LOOP%JUMP==0)
        {                              

		//pozycje---------------------------------------------------------------

		R_X.precision(nis_prec);
        for(int k=0;k<NATOMS;k++){R_X<<RX[k]<<" ";}
        R_X<<" "<<"\n";
        R_X.flush();	
		
		R_Y.precision(nis_prec);
        for(int k=0;k<NATOMS;k++){R_Y<<RY[k]<<" ";}
        R_Y<<" "<<"\n";
        R_Y.flush();
		
		R_Z.precision(nis_prec);
        for(int k=0;k<NATOMS;k++){R_Z<<RZ[k]<<" ";}
        R_Z<<" "<<"\n";
        R_Z.flush(); 
			 
        //----------------------------------------------------------------------
        cout.precision(5);
		}
//------------------------------------------------------------------------------
DENOM=DENOM+1;
}
//------------------------------------------------------------------------------        		
ostatnie_wartosci.precision(16);
ostatnie_wartosci<<VSH*158<<endl;
ostatnie_wartosci<<P0*42.133<<endl;
ostatnie_wartosci.flush();
//------------------------------------------------------------------------------
R_X.close(); R_Y.close(); R_Z.close();  przebiegi.close(); 		
//------------------------------------------------------------------------------

fstream plik_poloz(nazwa_ostatnie_pozycje, ios::out );
    if (plik_poloz.good()) 
	{  
            for(int i=0;i<NATOMS;i++){
                plik_poloz.precision(16);
                plik_poloz << RX[i]<< " \t " << RY[i] << " \t " << RZ[i] << " \t " << VX[i]<< " \t " << VY[i] << " \t " << VZ[i] << " \t " << RX_eq[i] << " \t " << RY_eq[i] << " \t " << RZ_eq[i] <<"\n";
                plik_poloz.flush();
            }

            plik_poloz.close();
    }
    
    fstream reszta_wejscie_new("reszta_wejscie_new.txt", ios::out );
    if (plik_poloz.good()) 
	{  
            
                reszta_wejscie_new.precision(16);
				reszta_wejscie_new <<VSH<<"\n";
				reszta_wejscie_new <<P0<<"\n";
				reszta_wejscie_new <<P0_K;
				reszta_wejscie_new <<VSH<<"\n";
				reszta_wejscie_new <<P0<<"\n";
				reszta_wejscie_new <<P0_K<<"\n";
				reszta_wejscie_new <<RUN<<"\n";
				reszta_wejscie_new <<cleaning<<"\n";
				reszta_wejscie_new <<a1<<"\n";
				reszta_wejscie_new <<a2<<"\n";
				reszta_wejscie_new <<BARO4<<"\n";
				reszta_wejscie_new <<BARO5<<"\n";
				reszta_wejscie_new <<beta<<"\n";
				reszta_wejscie_new <<READ<<"\n";
				reszta_wejscie_new <<wczytaj_sciany<<"\n";
				reszta_wejscie_new <<granica_odczytu_scian<<"\n";
				reszta_wejscie_new <<GAP<<"\n";
				reszta_wejscie_new.flush();
            reszta_wejscie_new.close();
    }
//------------------------------------------------------------------------------

                for (int i=0; i<RANG; i++) 
			{
				NWGR[i]=WGR[i];
				NGR[i]=GR[i];
				
                TT[i]=i*pow(bin,-1);
                
				NWGR[i]=NWGR[i]*pow(RUN,-1)*pow(NATOMS,-1);
				NGR[i]=NGR[i]*pow(RUN,-1)*pow(NATOMS-2*wskA,-1);
                NWGR[i]=NWGR[i]*pow((4*M_PI*(TT[i]*TT[i])*rho*(pow(bin,-1))),-1);
                NGR[i]=NGR[i]*pow((4*M_PI*(TT[i]*TT[i])*rho*(pow(bin,-1))),-1);
                
                NhistZ[i]=histZ[i]*pow(RUN,-1);
            }

GR[0]=0; WGR[0]=0; NGR[0]=0; NWGR[0]=0;

for(int i=0;i<TGAP;i++)
{
	if(tabDenomRTZ[i]>0){MEAN_RTZ2[i]=SUM_RTZ[i]*pow(tabDenomRTZ[i],-1);}
	else{MEAN_RTZ2[i]=0;}
	
	if(tabDenomRTZ_X[i]>0){MEAN_RTZ_X2[i]=SUM_RTZ_X[i]*pow(tabDenomRTZ_X[i],-1);}
	else{MEAN_RTZ_X2[i]=0;}
	
	if(tabDenomRTZ_Y[i]>0){MEAN_RTZ_Y2[i]=SUM_RTZ_Y[i]*pow(tabDenomRTZ_Y[i],-1);}
	else{MEAN_RTZ_Y2[i]=0;}
	
	if(tabDenomRTZ_Z[i]>0){MEAN_RTZ_Z2[i]=SUM_RTZ_Z[i]*pow(tabDenomRTZ_Z[i],-1);}
	else{MEAN_RTZ_Z2[i]=0;}
	
	if(tabDenomRVX[i]>0){MEAN_RTVX2[i]=SUM_RTVX[i]*pow(tabDenomRVX[i],-1);}
	else{MEAN_RTVX2[i]=0;}
	
	if(tabDenomRVY[i]>0){MEAN_RTVY2[i]=SUM_RTVY[i]*pow(tabDenomRVY[i],-1);}
	else{MEAN_RTVY2[i]=0;}
	
	if(tabDenomRVY_pec[i]>0){MEAN_RTVY_pec2[i]=SUM_RTVY_pec[i]*pow(tabDenomRVY_pec[i],-1);}
	else{MEAN_RTVY_pec2[i]=0;}
	
	if(tabDenomRVZ[i]>0){MEAN_RTVZ2[i]=SUM_RTVZ[i]*pow(tabDenomRVZ[i],-1);}
	else{MEAN_RTVZ2[i]=0;}
}


for(int i=0;i<TGAP;i++)
{
//------------------------------------------------------------------------------	
v_prof.precision(16);
v_prof<<
/*1*/(i-TGAP*0.5)*resol<<			"  "<<
/*2*/MEAN_RTVX[i]<<	"  "<</*3*/MEAN_RTVY[i]<<		"  "<</*4*/MEAN_RTVZ[i]<<"  "<<"\n";
v_prof.flush();
//------------------------------------------------------------------------------	
T_prof.precision(16);
T_prof<<
/*1*/(i-TGAP*0.5)*resol<<"  "<</*2*/MEAN_RTZ_2[i]<<"\n";
T_prof.flush();
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
for (int i=0;i<RANG; i++)
{
H_prof.precision(16);
H_prof<</*1*/(i-0.5*RANG)*pow(bin,-1)<<"   "<</*2*/histZ[i]*pow(RUN-cleaning,-1)<<"\n";
H_prof.flush();
}
//------------------------------------------------------------------------------   

for (int i=0;i<RANG2; i++)
{
hist_P_L.precision(16);
hist_P_L<</*4*/(i)*pow(binP,-1)+Pb<<"   "<</*5*/histP[i]<<"   "<</*6*/(i)*pow(binL,-1)+L0<<"   "<</*7*/histL[i]<<"\n";
hist_P_L.flush();
}
//------------------------------------------------------------------------------

v_prof.close(); T_prof.close(); 
H_prof.close(); 	ostatnie_wartosci.close();
hist_P_L.close();

//------------------------------------------------------------------------------

cout<<"koniec"<<endl;
getch();     
}



