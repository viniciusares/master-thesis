function [r1,c1]=n3_r1c1Regular(X0,Xcj,Ni,Nj,NJ,Sj,SJ,D,v)
% (j-=j,j+=J) antes e depois canto

[Rcj,dRcjdni]=n4nucleosR(X0,Xcj,Ni,Nj,NJ,Sj,SJ,v);

[~,~,w,dwdn]=n4nucleos(X0,Xcj,Ni,D,v);

% Deslocamento/Deflexao
wcj=w;

% Rotacao
dwcjdni=-dwdn; %dwdn=dwdr*drdn Gri=-Gr

%-------Avaliacao componentes r1------------------------
r1(1,1)=Rcj;        % r2c2 usa o mesmo nucleo
    
r1(2,1)=dRcjdni;
    
%-------Avaliacao componentes c1------------------------
c1(1,1)=wcj;        % r2c2 usa o mesmo nucleo
    
c1(2,1)=dwcjdni;
