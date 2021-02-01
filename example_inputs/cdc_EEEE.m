% CDCs Placa de Kirchhoff Zeros, toda engastada
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0.3; h=0.01;
D=E*h^3/(12*(1-v^2));

% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1     m1,     w2      m2,     w3      m3
cdcLado1=[1   0       0       0       0       0       0]; % E
cdcLado2=[1   0       0       0       0       0       0]; % E
cdcLado3=[1   0       0       0       0       0       0]; % E
cdcLado4=[1   0       0       0       0       0       0]; % E

CDC=GeraCDC(n,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
CTOS(:,5:6)=[...
    1 0
    1 0
    1 0
    1 0];

cargcel=1;
Pcels(:,4)=cargcel;