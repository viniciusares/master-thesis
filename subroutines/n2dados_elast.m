% Dados Placa de Kirchhoff (Elasticidade)
    % elementos descontinuos

D=4; v=0.3;

%   num x   y   canto
GEO=[
    1   0   0   1
    2   3   0   0
    3   6   0   1
    4   6   3   0
    5   6   6   1
    6   3   6   0
    7   0   6   1
    8   0   3   0];

ELEMGEO=[
    1 1 2 3
    2 3 4 5
    3 5 6 7
    4 7 8 1];

% Matriz CANTOS, colunas: 1num 2geo 3antes 4depois
J=0; % inicializa contador de cantos
for I=1:size(GEO,1)         % percorre pontos geometricos
    if GEO(I,4)==1              % se geo for um canto:
        J=J+1;                  % ->conta esse canto
        CANTOS(J,1:2)=[J,I];    % ->e armazena na matriz cantos
        for K=1:size(ELEMGEO,1)     % percorre os elementos     
            if ELEMGEO(K,4)==CANTOS(J,2)    % se geo3 for cantoJ:
                CANTOS(J,3)=K;              % -> armazena como "antes"
            elseif ELEMGEO(K,2)==CANTOS(J,2)% se geo1 for cantoI:
                CANTOS(J,4)=K;              % -> armazena como "depois"
            end
        end
    end
end

for I=1:size(ELEMGEO,1) % percorre os elementos
    k=ELEMGEO(I,2:4);   % armazenas os 3 geos
    x1=GEO(k(1),2); y1=GEO(k(1),3); % extrai as coordenadas
    x2=GEO(k(2),2); y2=GEO(k(2),3);
    x3=GEO(k(3),2); y3=GEO(k(3),3);
    NOS(3*I-2:3*I,2:3)=n3nosC(x1,y1,x2,y2,x3,y3); % escreve nos de 3 em 3
end

% num, tipo, w1 m1, w2 m2, w3 m3
CDC=[1 2    300 0   300 0   300 0
    2 2     300 0   300 0   300 0
    3 2     350 0   400 0   350 0
    4 2     300 0   300 0   300 0];
% Tipos de CDC (codes)----Dados---------Incognitas--
% tipo 1: engastada     w=0 dwdn=0      Vn=? Mn=?
% tipo 2: apoiada       w=0 Mn=0        Vn=? dwdn=?
% tipo 3: guia          Vn=0 dwdn=0     w=? Mn=?
% tipo 4: livre         Vn=0 Mn=0       w=? dwdn=?

% CDC cantos col5:1-preso 2-solto, col6:valorCDC
% dado 1-preso: wc || dado 2-solto: Rc
CANTOS(:,5:6)=[...
    1 0
    1 0
    1 0
    1 0];

P_int=[1 1 3
        2 1.5 3
        3 2 3
        4 3 3
        5 4 3
        6 4.5 3
        7 5 3
        8 3 1
        9 3 1.5
        10 3 2
        11 3 3
        12 3 4
        13 3 4.5
        14 3 5];