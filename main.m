clear all
close all
clc
format long

% Versao Dezembro (09/01/2016)
% Programa MEC Placa - Vinicius Emanoel Ares 2015
% elementos quadraticos, descontinuos, retos

% DADOS DE ENTRADA (Pre-Processamento) %=============================
% nel=1; Lado=1; % num. de Elem./Lado, dados placa quadrada
nelX=2;nelY=2;LadoX=1;LadoY=1; % Dados placa retangular
cruzougrade=2; % P_ints: 1-Cruz_EALA, 2-Grade, 3-Cruz_AAAA
%------------------------------------------------
%geometria % N Elem./Lado
geometriaPlacaViga

%------------------------------------------------
cdc_LALE
% cdc_LELE
% cdc_ALLA
% cdc_LALA
% cdc_AAAAxy
% cdc_ELLL       % com_cargcel, Eixo da Viga: Y
% cdc_LLLE       % com_cargcel, Eixo da Viga: X
% cdc_AAAA       % carg_uniforme
% cdc_EALA       % Engast-Apoiado-Livre-Apoiado
%cdc_AALA       % Apoiado-Apoiado-Livre-Apoiado
%cdc_EEEE       % com_cargcel
%cdc_ALAL100    % placa_inclinada
%cdc_EEEE100    % SO SERVE P 1E/L e 6x6m
%cdc_EEEE100x3  % SO SERVE P 3E/L e 6x6m
%------------------------------------------------

% PRINCIPAL (Solver) %===============================================
[H1,G1]=n2montaH1G1sing5A(NOS,ELGEO,GEO,D,v);
[R1,C1]=n2montaR1C1(NOS,ELGEO,GEO,CTOS,D,v);
[H2,G2]=n2montaH2G2sing3(CTOS,ELGEO,GEO,D,v);
%[H2,G2]=n2montaH2G2sing6(CTOS,ELGEO,GEO,D,v);
[R2,C2]=n2montaR2C2(CTOS,ELGEO,GEO,D,v,H2);

% Junta as submatrizes para montar as matrizes globais
H=[ H1 R1
    H2 R2];
G=[ G1 C1
    G2 C2];

[A,b]=n2aplica_CDC(CDC,H,G,CTOS); % Aplica as CDCs
Cc=CargasDominio(NOS,ELGEO,Pcels,Celulas,D,2); % tipo 2 por no
Cbn=CargasDominio(GEO(CTOS(:,2),1:3),0,Pcels,Celulas,D,1);
C=[Cc;Cbn];
x=A\(b+C); % Calcula as variaveis desconhecidas
%xD=VectorX_Daros;


% ORGANIZA E ESCREVE SOLUCAO
format short
[w,Vn,dwdn,Mn,wc,Rc]=n2monta_respostas(x,CDC,CTOS); % corrigido 1->I
w_dwdn=[w' dwdn'] %#ok
Vn_Mn=[Vn' Mn'] %#ok
wc=wc' %#ok
Rc=Rc' %#ok

%======================================================================
% POS-PROCESSAMENTO
if cruzougrade==2 % 2-Grade
    clear P_int % Apaga P_int CRUZ e define Grade
    xmin=0.10*LadoX; ymin=0.10*LadoY;   %CUIDADO!!!
    xmax=0.90*LadoX; ymax=0.90*LadoY;   %CUIDADO!!!
    xnpts=9; ynpts=9;                   %CUIDADO!!!
    P_int=FgeraGrade(xmin,ymin,xmax,ymax,xnpts,ynpts);
end

% CALCULA E ESCREVE OS VALORES INTERNOS
format short
[H1i,G1i,R1i,C1i,u,q]=wIntern(P_int,ELGEO,GEO,CTOS,w,dwdn,Vn,Mn,wc,Rc,D,v);
Cint=CargasDominio(P_int,0,Pcels,Celulas,D,1);
w_int=-[H1i R1i]*u'+[G1i C1i]*q'+Cint;

%compara=[VectorP_Daros , C]
%(VectorP_Daros./C-ones(28,1))*100

variosgraficos
Iycte=ynpts/2+0.5; % y=0.5 Meio-caminho entre E e L
deflsFx=Wgrid(Iycte,:);
format short
% xgrid(10)=1;
% deflsFx(10)=w_dwdn(5,1);
(deflsFx'-LALEdFxDarosSimples1)./LALEdFxDarosSimples1*100
figure(5)
plot(xgrid,LALEdFxDarosSimples1,'b')    % daros simples
hold on
plot(xgrid,deflsFx,'b--')               % prog simples
plot(xgrid,LALEdFxDarosStretch1,'r')    % daros stretch
plot(xgrid,LALEdFxStretch1,'r--')       % prog stretch
legend('DarosSimples','ProgSimples','DarosStretch','ProgStretch')