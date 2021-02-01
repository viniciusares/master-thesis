% CDCs Placa de Kirchhoff Engasrigida 0-200mm x1 (1elem/lado)
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
cdcLado1=[1   0.000   -teta   0.000   -teta   0.000   -teta]; % E
cdcLado2=[1   0.3/9   0       0.9/9   0       1.5/9   0    ]; % E
cdcLado3=[1   0.200   teta    0.200   teta    0.200   teta ]; % E
cdcLado4=[1   1.5/9   0       0.9/9   0       0.3/9   0    ]; % E

n=1; % numero de elem./lado (NAO MUDAR NESSE CASO)
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