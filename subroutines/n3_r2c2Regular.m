function [r2,c2]=n3_r2c2Regular(X0,Xcj,Nj,NJ,Sj,SJ,D,v)
% (j-=j,j+=J) antes e depois canto

Ni=[1;1]; % nao usado

[Rcj,~]=n4nucleosR(X0,Xcj,Ni,Nj,NJ,Sj,SJ,v);

[~,~,w,~]=n4nucleos(X0,Xcj,Ni,D,v);

wcj=w;

%-------Avaliacao componentes r2------------------------
r2=Rcj; % Forca de Canto %igual a n3_r1c1Regular
    
%-------Avaliacao componentes c2------------------------
c2=wcj; % Deslocamento/Deflexao %igual a n3_r1c1Regular