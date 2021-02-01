function GEO=GeraGEOxy(nelX,nelY,Lx,Ly)
% Essa rotina gera a matriz GEO apenas para o caso de placa (Lx por Ly)
% nelX,nelY eh o numero de elementos (quadraticos) por lado que se deseja

denX=2*nelX; denY=2*nelY; % denominadores


lado1=zeros(2*nelX,4); %num x y canto
for I=1:2*nelX
    lado1(I,1)=I;               % numero do ponto
    lado1(I,2)=Lx*(I-1)/denX;   % coordenada x (crescente)
    lado1(I,3)=0;               % coordenada y=0
end
lado1(1,4)=1; % o primeiro ponto eh Canto

lado2=zeros(2*nelY,4);
for I=1:2*nelY
    lado2(I,1)= 2*nelX + I;     % numero do ponto
    lado2(I,2)=Lx;              % coordenada x=L
    lado2(I,3)=Ly*(I-1)/denY;   % coordenada y (crescente)
end
lado2(1,4)=1; % o primeiro ponto eh Canto

lado3=zeros(2*nelX,4);
for I=1:2*nelX
    lado3(I,1)= 2*nelX + 2*nelY + I;    % numero do ponto
    lado3(I,2)=Lx*(2*nelX+1-I)/denX;    % coordenada x (decrescente)
    lado3(I,3)=Ly;                      % coordenada y=L
end
lado3(1,4)=1; % o primeiro ponto eh Canto

lado4=zeros(2*nelY,4);
for I=1:2*nelY
    lado4(I,1)= 4*nelX + 2*nelY + I;    % numero do ponto
    lado4(I,2)=0;                       % coordenada x=0
    lado4(I,3)=Ly*(2*nelY+1-I)/denY;    % coordenada y (decrescente)
end
lado4(1,4)=1; % o primeiro ponto eh Canto

GEO=[lado1
    lado2
    lado3
    lado4];