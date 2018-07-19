//geometria ukladu--------------------------------------------------------------

mass=1;
KB=1;
g=3*NATOMS-3;
DIMRATE=NUNITZ*pow(NUNITX,-1);

rho=1.1;
READ=1;//1-dane, 2-dane scian, 3-blok
wczytaj_sciany=0;

LEY=4;

ZerowaniePedu=1;//dot. READ=1
ZerowanieSrodkaM=1;//dot. READ=1
zerowanieSM2=0;

XPERIOD=1;
YPERIOD=1;
ZPERIOD=0;

GAP=0.5;
CENTGAP=0;

//szczegoly symulacji-----------------------------------------------------------

RUN=12000;
DT=0.0025;
cleaning=2000;
n_0=10; //odstep miedzy zmiana warunkow symulacji
JUMP=200; //dawac jak najwieksze - pozycje i chwilowe histogramy
JUMP2=200; //przebiegi czasowe
t_test=3000;
odczyt_wejscie=1;//odczyt,
dataLines=10;
calcGR=0;
StartVelocities=1;
CONSTR=1; //0 - uwiazanie atomow scian w danych punktach, 1 - sciany ruchome, 2 - sciany wolne
DIFINT=1; //rozne oddzialywania miedzy grupami atomow 
DIFCON=1; //zastosowanie roznych warunkow dla roznych warstw (cisnienie, scinanie, temperatura)
faktor_k=0.001;
N1=12; N2=6; 
k4=faktor_k*5000; 
k6=faktor_k*5000000; //do potencjalu anharmonicznego PH
potencjal=0; //0 - LJ, 1 - WCA;
RCUT=2.5;
WCA=pow(2,pow(6,-1));//odciecie dla WCA
SKIN=1.5;
R3=1.5;
KR3=0.0001;
WR3=12;

TETF=0; //TETF - potencjal PH
TETF2=0; //TETF2 - potecjal 0.5*k*r^2
PAIR_TET=0;
//temperatura-------------------------------------------------------------------

TERM=1;
T0=1;
T0_up=T0;
T0_down=T0;
scale_vel_w=2; //1 - skalowanie predkosci w objetosci, 2 - skalowanie predkosci w scianach, 0-nic

korekta_pedu_scian=1;

freq_scale=1;
zmiana_warunkow_T=0; //zmiana z grzania jednorodnego na w scianach w chwili n_0

termWall=1;//0 - grzanie jednorodne, 1 - wspolna temperatura dla obu scian, 2 - oddzielne temperatury dla obu scian
termX=1;
termY=1;
termZ=1;

//scinanie----------------------------------------------------------------------

SHEAR=1;
automat_sliding=1;
sliding=1; //rozpedzanie ukladu
deltaV=1;
VSH=0; //poczatkowa predkosc scinania
VSH_K=0;

SH1=1;
SH2=0;
SH3=0;
FY0=0.5;

//cisnienie---------------------------------------------------------------------

BAR=0;
BAR_CHOICE=0;//1-BH; 2-DHO Verlet; 3-DHO RK4; 4-HO; 0-nic
automat_pressing=1;
pressing=0; //1-sciskanie ukladu; 2-rozprezanie;
cisnienieNaturalne=0;

deltaP=1; //zawsze dodatnie! (sciskanie pressing 1, rozprezanie pressing 2)
deltaJumpP=-2;

BARO4=1;
BARO5=1;

sygnal=4; //1-Pt; 2-VAP; 3-PVIR1; 4-PW;
P0_K=6;
P0=6;

Qp=1;
omega=2;
beta=1;
gam1=0.5*P0*BX*BX*VSH*(1/(VSH*VSH+1));
gam2=250;
//gam1=0;

jump_P0=0; //skok docelowego cisnienia

mop_calc=0;
PN_calc=0;
tauH=pow(omega,-1);

//histogramy--------------------------------------------------------------------
L0=NUNITZ; //oczekiwana srednia szerokosc ukladu,

binL=300; //histogram L

Pb=0.75*P0_K;
ZP=5;//okresla okresla zasieg liczenia cisnienia
binT=100;
binD=100;
binP=100; //histogram P

bin=10; // szerokosc binow HR i GR
binRZ=1;
binHV=25;

resol=1;

