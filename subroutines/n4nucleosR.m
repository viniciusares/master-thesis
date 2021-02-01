function [Rcj,dRcjdni]=n4nucleosR(X0,X,Ni,Nj,NJ,Sj,SJ,v)

% Parametros
x0=X0(1); y0=X0(2); x=X(1); y=X(2);
Rx=x-x0; Ry=y-y0; % R deve ser diferente de 0, nao singular
R=sqrt(Rx^2+Ry^2); % distancia entre o ponto fonte e ponto campo
Gr=[Rx/R,Ry/R]; % G=gradiente [ r,x , r,y ] (Gr*N=cosb=drdn)

% Forca de Canto  % usada em n3_r1c1Regular e n3_r2c2Regular
Rcj=-(1-v)/(4*pi)*((Gr*NJ)*(Gr*SJ)-(Gr*Nj)*(Gr*Sj));

% Derivada da Forca de Canto  % usada em n3_r1c1Regular
dRcjdni=-(1-v)/(4*pi*R)*(...
+( (Gr*SJ)*((Gr*NJ)*(Gr*Ni)-(Ni'*NJ))...
  +(Gr*NJ)*((Gr*SJ)*(Gr*Ni)-(Ni'*SJ)) )...
-( (Gr*Sj)*((Gr*Nj)*(Gr*Ni)-(Ni'*Nj))...
  +(Gr*Nj)*((Gr*Sj)*(Gr*Ni)-(Ni'*Sj)) )...
);

