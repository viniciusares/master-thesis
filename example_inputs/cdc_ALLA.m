% CDCs Placa de Kirchhoff Caso LLAA carga uniforme
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0; h=0.010;
D=E*h^3/(12*(1-v^2));

% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1 m1, w2 m2, w3 m3
cdcLado1=[2   0       0       0       0       0       0]; % A y=0
cdcLado2=[4   0       0       0       0       0       0]; % L x=Lx
cdcLado3=[4   0       0       0       0       0       0]; % L y=Ly
cdcLado4=[2   0       0       0       0       0       0]; % A x=0

CDC=GeraCDCxy(nelX,nelY,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
CTOS(:,5:6)=[...
    1 0
    1 0
    2 0
    1 0];

cargcel=1;
Pcels(:,4)=cargcel;