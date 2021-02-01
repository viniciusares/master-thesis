% Dados Placa de Kirchhoff (Elasticidade) Engastado-Apoiado
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
cdcLado4=[1   0       0       0       0       0       0]; % E x=0

CDC=GeraCDCxy(nelX,nelY,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
cargpontual=0; % Forca Pontual
CTOS(:,5:6)=[...
    1 0
    1 0
    1 0
    1 0];

cargcel=1; % Forca de Superficie
Pcels(:,4)=cargcel;

LALEdFxDarosSimples1=1e-7*[
    0.315
    1.04
    1.89
    2.64
    3.125
    3.24
    2.94
    2.24
    1.215];

LALEdFxStretch1=1.0e-06*[
   0.027523835701137
   0.088287622061122
   0.156677848934666
   0.214555790894244
   0.250079819947064
   0.256387785574771
   0.231005972275853
   0.175333533856562
   0.094409114162157];

LALEdFxDarosStretch1=1e-7*[
    0.259546
    0.840374
    1.5029
    2.07278
    2.43036
    2.50351
    2.26323
    1.72195
    0.934201];