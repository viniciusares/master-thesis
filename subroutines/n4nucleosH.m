function [dVndni,dMndni,dwdni,d2wdndni]=n4nucleosH(X0,X,Ni,N,D,v)

% Parametros
x0=X0(1); y0=X0(2); x=X(1); y=X(2);
Rx=x-x0; Ry=y-y0; % R deve ser diferente de 0, nao singular
R=sqrt(Rx^2+Ry^2); % distancia entre o ponto fonte e ponto campo
Gr=[Rx/R,Ry/R]; % G=gradiente [ r,x , r,y ] (Gr*N=cosb=drdn)

% Derivada da Forca Cortante em relacao aa fonte %ok29/04/15
dVndni=...
-(Gr*N)*(Gr*Ni)/(4*pi*R^2)*((4+2*(1-v)*(2*(Gr*N)^2-1))+4*(1-v)*(Gr*N)^2)...
        +(N'*Ni)/(4*pi*R^2)*((2+(1-v)*(2*(Gr*N)^2-1))+4*(1-v)*(Gr*N)^2);
%OK com Nucleos da Sol. Fund. de Placas 01/03/15
    
% Derivada do Momento %ok29/04/15
dMndni=1/(4*pi*R)*((1+v)*(Gr*Ni)-2*(1-v)*((Gr*N)^2*(Gr*Ni)-(Gr*N)*(N'*Ni)));
%OK com Nucleos da Sol. Fund. de Placas 01/03/15
% Difere de Nucleos da Sol. Fund. usando Artigo de Stern 27/02/15 (1+v)
%dMndni=-dMndni; %Divergencia com Mathematica Daros 11/05/15

% Rotacao em relacao ao no fonte
dwdr=1/(8*pi*D)*(2*R*log(R)+R);
dwdni=-dwdr*(Gr*Ni); %dwdn=dwdr*drdn %24/04/15 Gri=-Gr

% Derivada da Rotacao %ok29/04/15
d2wdndni=-1/(8*pi*D)*((2*log(R)+1)*(N'*Ni)+2*(Gr*N)*(Gr*Ni));
%OK com Nucleos da Sol. Fund. usando Artigo de Stern 27/02/15
%OK com Nucleos da Sol. Fund. de Placas 01/03/15