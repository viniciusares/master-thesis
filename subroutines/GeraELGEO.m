function ELGEO = GeraELGEO(nel)

% col1: num do elemento, cols 2 a 4: pontos GEO 1 a 3
ELGEO=zeros(4*nel,4);

for I=1:4*nel
    ELGEO(I,1)=I;       % col1: num do elemento
    ELGEO(I,2)=2*I-1;   % col2: ponto GEO 1
    ELGEO(I,3)=2*I;     % col3: ponto GEO 2
    ELGEO(I,4)=2*I+1;   % col4: ponto GEO 3
end

ELGEO(4*nel,4)=1; % sobrescreve ultimo ponto para eliminar ponto que
% excederia o numero de pontos GEO
% Sobrescreve 1, dessa forma fechando o "quadrado" da placa.

% Exemplo de como fica ELGEO apos pronto:
% ELGEO=[
%     1 1 2 3
%     2 3 4 5
%     3 5 6 7
%     4 7 8 1];