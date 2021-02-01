% Dados Placa de Kirchhoff (Elasticidade) Engastado-Livre
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0; h=0.010;
D=E*h^3/(12*(1-v^2));

carg=0; % Carga Linear
% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1 m1, w2 m2, w3 m3
cdcLado1=[1   0       0       0       0       0       0]; % E y=0
cdcLado2=[4   0       0       0       0       0       0]; % L x=Lx
cdcLado3=[4   carg    0       carg    0       carg    0]; % L y=Ly
cdcLado4=[4   0       0       0       0       0       0]; % L x=0

CDC=GeraCDCxy(nelX,nelY,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
cargpontual=0;
CTOS(:,5:6)=[...
    1 0
    1 0
    2 cargpontual
    2 cargpontual];

cargcel=1; % Carga de Superficie
Pcels(:,4)=cargcel;