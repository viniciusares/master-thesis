% geometria placa-viga 10x1

% Colunas GEO: num   x   y   canto
GEO=GeraGEOxy(nelX,nelY,LadoX,LadoY); % nel = numero Elem./Lado

% num x   y (col4:cargas dadas em cdc)
Pcels=[
    1 0             0
    2 0.5*LadoX     0
    3 0.5*LadoX     0.5*LadoY
    4 0             0.5*LadoY
    5 1*LadoX       0
    6 1*LadoX       0.5*LadoY
    7 1*LadoX       1*LadoY
    8 0.5*LadoX     1*LadoY
    9 0             1*LadoY];

ELGEO=GeraELGEO((nelX+nelY)/2); % n = numero Elem./Lado

% Incidencia celular, colunas: n_cel + 4 p_cels
Celulas=[
    1 1 2 3 4
    2 2 5 6 3
    3 3 6 7 8
    4 4 3 8 9];

% Matriz CTOS(cantos), colunas: 1num 2geo 3elemantes 4elemdepois
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
    1 0.1*LadoX 0.5*LadoY
    2 0.2*LadoX 0.5*LadoY
    3 0.3*LadoX 0.5*LadoY
    4 0.4*LadoX 0.5*LadoY
    5 0.5*LadoX 0.5*LadoY
    6 0.6*LadoX 0.5*LadoY
    7 0.7*LadoX 0.5*LadoY
    8 0.8*LadoX 0.5*LadoY
    9 0.9*LadoX 0.5*LadoY
    10 0.5*LadoX 0.1*LadoY
    11 0.5*LadoX 0.2*LadoY
    12 0.5*LadoX 0.3*LadoY
    13 0.5*LadoX 0.4*LadoY
    14 0.5*LadoX 0.5*LadoY
    15 0.5*LadoX 0.6*LadoY
    16 0.5*LadoX 0.7*LadoY
    17 0.5*LadoX 0.8*LadoY
    18 0.5*LadoX 0.9*LadoY];