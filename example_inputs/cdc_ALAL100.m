% CDCs Placa de Kirchhoff 0-200mm
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0.3; h=0.006;
D=E*h^3/(12*(1-v^2));

% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1     m1,     w2      m2,     w3      m3
cdcLado1=[2   0.000   0       0.000   0       0.000   0]; % A
cdcLado2=[4   0       0       0       0       0       0]; % L
cdcLado3=[2   0.200   0       0.200   0       0.200   0]; % A
cdcLado4=[4   0       0       0       0       0       0]; % L

CDC=GeraCDC(n,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
CTOS(:,5:6)=[...
    1 0
    1 0
    1 0.200
    1 0.200];

cargcel=0;
Pcels(:,4)=cargcel;