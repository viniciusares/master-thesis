% CDCs Placa de Kirchhoff Caso LALA carga uniforme
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0; h=0.010;
D=E*h^3/(12*(1-v^2));

% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1 m1, w2 m2, w3 m3
cdcLado1=[4   0       0       0       0       0       0]; % L y=0
cdcLado2=[2   0       0       0       0       0       0]; % A x=Lx
cdcLado3=[4   0       0       0       0       0       0]; % L y=Ly
cdcLado4=[2   0       0       0       0       0       0]; % A x=0

CDC=GeraCDCxy(nelX,nelY,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
CTOS(:,5:6)=[...
    1 0
    1 0
    1 0
    1 0];

cargcel=1;
Pcels(:,4)=cargcel;

LALAdFxStretch1=1.0e-06*[
   0.152664632444412
   0.289919281392300
   0.395987076075741
   0.462916277296637
   0.485780050945541
   0.462916277296636
   0.395987076075741
   0.289919281392300
   0.152664632444412];

LALAdeflsFxDarosStretch1=[
    1.5345e-7
    2.89519e-7
    3.95307e-7
    4.62117e-7
    4.84933e-7
    4.62117e-7
    3.95307e-7
    2.89519e-7
    1.5345e-7];

LALAdFxDarosSimples1=1e-7*[
    2.4525
    4.64
    6.3525
    7.44
    7.8125
    7.44
    6.3525
    4.64
    2.4525];