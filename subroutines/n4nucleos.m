function [Vn,Mn,w,dwdn]=n4nucleos(X0,X,N,D,v)

% Parametros
x0=X0(1); y0=X0(2); x=X(1); y=X(2);
Rx=x-x0; Ry=y-y0; % R deve ser diferente de 0, nao singular
R=sqrt(Rx^2+Ry^2); % distancia entre o ponto fonte e ponto campo
Gr=[Rx/R,Ry/R]; % G=gradiente [ r,x , r,y ] (Gr*N=cosb=drdn)

% Forca Cortante
Vn=-(Gr*N)/(4*pi*R)*(2+(1-v)*(2*(Gr*N)^2-1)); %ok29/04/15
%OK com Nucleos da Sol. Fund. usando Artigo de Stern 27/02/15
% difere apenas no raio de curvatura
% Difere de Nucleos da Sol. Fund. de Placas 01/03/15 (1-v)

% Momento (fletor)
Mn=-(1+v)/(4*pi)*(1+log(R))-(1-v)/(8*pi)*(2*(Gr*N)^2-1); %ok29/04/15
%OK com Nucleos da Sol. Fund. usando Artigo de Stern 27/02/15
%OK com Nucleos da Sol. Fund. de Placas 01/03/15
%Mn=-Mn; %Divergencia com Mathematica Daros 11/05/15

% Deslocamento/Deflexao
w=1/(8*pi*D)*R^2*log(R); %ok29/04/15
%OK com Nucleos da Sol. Fund. usando Artigo de Stern 27/02/15
%OK com Nucleos da Sol. Fund. de Placas 01/03/15

% Rotacao normal do bordo
dwdr=1/(8*pi*D)*(2*R*log(R)+R);
dwdn=dwdr*(Gr*N); %dwdn=dwdr*drdn %ok29/04/15
%OK com Nucleos da Sol. Fund. usando Artigo de Stern 27/02/15
%OK com Nucleos da Sol. Fund. de Placas 01/03/15