% CDCs Placa de Kirchhoff Caso LELE carga uniforme
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
cdcLado2=[1   0       0       0       0       0       0]; % E x=Lx
cdcLado3=[4   0       0       0       0       0       0]; % L y=Ly
cdcLado4=[1   0       0       0       0       0       0]; % E x=0

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

LELEdFxDaros10=[
    0.000202500000000
   0.000640000000000
   0.001102500000000
   0.001440000000000
   0.001562500000000
   0.001440000000000
   0.001102500000000
   0.000640000000000
   0.000202500000000];

LELEdFxDarosSimples1=[
    2.025e-8
    6.4e-8
    1.1025e-7
    1.44e-7
    1.5625e-7
    1.44e-7
    1.1025e-7
    6.4e-8
    2.025e-8];

LELEdFxStretch1=1.0e-06*[
   0.019134989049715
   0.059423990722152
   0.101088544618623
   0.130987126750403
   0.141762891717713
   0.130987126750402
   0.101088544618624
   0.059423990722152
   0.019134989049715];

LELEdFxDarosStretch1=[
    1.81604e-8
    5.66314e-8
    9.66276e-8
    1.25486e-7
    1.35901e-7
    1.25486e-7
    9.66276e-8
    5.66314e-8
    1.81604e-8];