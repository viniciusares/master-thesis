% CDCs Placa de Kirchhoff Engasrigida 0-200mm x3 (3elem/lado)
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0.3; h=0.006;
D=E*h^3/(12*(1-v^2)); %Pa.m3=N.m

teta=0.200/6;
% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1     dwdn1,  w2      dwdn2,  w3      dwdn3
CDC=[1 1    0.000   -teta   0.000   -teta   0.000   -teta
    2 1     0.000   -teta   0.000   -teta   0.000   -teta
    3 1     0.000   -teta   0.000   -teta   0.000   -teta
    4 1     0.1/9   0       0.3/9   0       0.5/9   0
    5 1     0.7/9   0       0.9/9   0       1.1/9   0
    6 1     1.3/9   0       1.5/9   0       1.7/9   0
    7 1     0.200   teta    0.200   teta    0.200   teta
    8 1     0.200   teta    0.200   teta    0.200   teta
    9 1     0.200   teta    0.200   teta    0.200   teta
    10 1    1.7/9   0       1.5/9   0       1.3/9   0
    11 1    1.1/9   0       0.9/9   0       0.7/9   0
    12 1    0.5/9   0       0.3/9   0       0.1/9   0];

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
CTOS(:,5:6)=[...
    1 0
    1 0
    1 0.200
    1 0.200];

cargcel=0;
Pcels(:,4)=cargcel;