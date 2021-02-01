% Dados Placa de Kirchhoff (Elasticidade) Engastado-Livre
% Define: v, D, CDC; Complementa: CTOS, Pcels
E=200e9; v=0; h=0.010;
D=E*h^3/(12*(1-v^2));

carg=0; % Forca Linear
% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?
% num, tipo, w1 m1, w2 m2, w3 m3
cdcLado1=[4   0       0       0       0       0       0]; % L y=0
cdcLado2=[4   carg    0       carg    0       carg    0]; % L x=Lx
cdcLado3=[4   0       0       0       0       0       0]; % L y=Ly
cdcLado4=[1   0       0       0       0       0       0]; % E x=0

CDC=GeraCDCxy(nelX,nelY,cdcLado1,cdcLado2,cdcLado3,cdcLado4);

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado preso: wc || dado solto: Rc
cargpontual=0; % Forca Pontual
CTOS(:,5:6)=[...
    1 0
    2 cargpontual
    2 cargpontual
    1 0];

cargcel=1; % Forca de Superficie
Pcels(:,4)=cargcel;

LLLEdFxDarosSimples1=[
    1.4025e-7
    5.24e-7
    1.10025e-6
    1.824e-6
    2.65625e-6
    3.564e-6
    4.52025e-6
    5.504e-6
    6.50025e-6
7.5e-6];

% LLLEdFxStretch1=1.0e-05*[
%    0.030668897778103
%    0.067945877380817
%    0.105344104023865
%    0.136978323944080
%    0.157667360813906
%    0.162993846366436
%    0.149496459999400
%    0.115337329276184
%    0.061437233203333
% -5.084e-9];

LLLEdFxDarosStretch1=1e-6*[
    0.12204
    0.451871
    0.940447
    1.54568
    2.23218
    2.97105
    3.73972
    4.52179
    5.30695
    6.0909];