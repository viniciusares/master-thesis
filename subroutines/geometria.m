% Geometria Placa de Kirchhoff, elementos descontinuos
% define: GEO, Pcels, ELGEO, Celulas, CTOS, NOS, P_int

% Colunas GEO: num   x   y   canto
GEO=GeraGEO(nel,Lado); % n = numero Elem./Lado

% num x   y (col4:cargas dadas em cdc)
Pcels=[
    1 0         0
    2 0.5*Lado  0
    3 0.5*Lado  0.5*Lado
    4 0         0.5*Lado
    5 1*Lado    0
    6 1*Lado    0.5*Lado
    7 1*Lado    1*Lado
    8 0.5*Lado  1*Lado
    9 0         1*Lado];

ELGEO=GeraELGEO(nel); % n = numero Elem./Lado

% Incidencia celular, colunas: n_cel + 4 p_cels
Celulas=[
    1 1 2 3 4
    2 2 5 6 3
    3 3 6 7 8
    4 4 3 8 9];

% Matriz CTOS(cantos), colunas: 1num 2geo 3antes 4depois
J=0; % inicializa contador de cantos
for I=1:size(GEO,1)         % percorre pontos geometricos
    if GEO(I,4)==1              % se geo for um canto:
        J=J+1;                  % ->conta esse canto
        CTOS(J,1:2)=[J,I];    % ->e armazena na matriz cantos
        for K=1:size(ELGEO,1)     % percorre os elementos     
            if ELGEO(K,4)==CTOS(J,2)    % se geo3 for cantoJ:
                CTOS(J,3)=K;              % -> armazena elem como "antes"
            elseif ELGEO(K,2)==CTOS(J,2)% se geo1 for cantoJ:
                CTOS(J,4)=K;              % -> armazena elem como "depois"
            end
        end
    end
end

% Gerador da matriz NOS
NOS=zeros(3*size(ELGEO,1),3);
for I=1:size(ELGEO,1) % percorre os elementos
    k=ELGEO(I,2:4);   % armazena os 3 geos
    x1=GEO(k(1),2); y1=GEO(k(1),3); % extrai as coordenadas
    x2=GEO(k(2),2); y2=GEO(k(2),3);
    x3=GEO(k(3),2); y3=GEO(k(3),3);
    NOS(3*I-2,1)=3*I-2; NOS(3*I-1,1)=3*I-1; NOS(3*I,1)=3*I;
    NOS(3*I-2:3*I,2:3)=n3nosC(x1,y1,x2,y2,x3,y3); % escreve nos de 3 em 3
end

P_int=[
    1 0.1*Lado 0.5*Lado
    2 0.2*Lado 0.5*Lado
    3 0.3*Lado 0.5*Lado
    4 0.4*Lado 0.5*Lado
    5 0.5*Lado 0.5*Lado
    6 0.6*Lado 0.5*Lado
    7 0.7*Lado 0.5*Lado
    8 0.8*Lado 0.5*Lado
    9 0.9*Lado 0.5*Lado
    10 0.5*Lado 0.1*Lado
    11 0.5*Lado 0.2*Lado
    12 0.5*Lado 0.3*Lado
    13 0.5*Lado 0.4*Lado
    14 0.5*Lado 0.5*Lado
    15 0.5*Lado 0.6*Lado
    16 0.5*Lado 0.7*Lado
    17 0.5*Lado 0.8*Lado
    18 0.5*Lado 0.9*Lado];