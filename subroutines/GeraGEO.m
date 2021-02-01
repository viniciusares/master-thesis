function GEO=GeraGEO(n,L)
% Essa rotina gera a matriz GEO apenas para o caso de placa LxL
% n eh o numero de elementos (quadraticos) por lado que se deseja

den=2*n; % denominador

lado1=zeros(2*n,4); %num x y canto
for I=1:2*n
    lado1(I,1)=I;           % numero do ponto
    lado1(I,2)=L*(I-1)/den;   % coordenada x (crescente)
    lado1(I,3)=0;           % coordenada y=0
end
lado1(1,4)=1; % o primeiro ponto eh Canto

lado2=zeros(2*n,4);
for I=1:2*n
    lado2(I,1)=2*n+I;       % numero do ponto
    lado2(I,2)=L;           % coordenada x=L
    lado2(I,3)=L*(I-1)/den;   % coordenada y (crescente)
end
lado2(1,4)=1; % o primeiro ponto eh Canto

lado3=zeros(2*n,4);
for I=1:2*n
    lado3(I,1)=4*n+I;           % numero do ponto
    lado3(I,2)=L*(2*n+1-I)/den;   % coordenada x (decrescente)
    lado3(I,3)=L;               % coordenada y=L
end
lado3(1,4)=1; % o primeiro ponto eh Canto

lado4=zeros(2*n,4);
for I=1:2*n
    lado4(I,1)=6*n+I;           % numero do ponto
    lado4(I,2)=0;               % coordenada x=0
    lado4(I,3)=L*(2*n+1-I)/den;   % coordenada y (decrescente)
end
lado4(1,4)=1; % o primeiro ponto eh Canto

GEO=[lado1
    lado2
    lado3
    lado4];