function [H1,G1]=n2montaH1G1sing5A(NOS,ELGEO,GEO,D,v)

n_el = size(ELGEO,1); % Numero total de elementos
n_nos = size(NOS,1); % Numero total de nos

% Inicializacao das Matrizes
H1 = zeros(2*n_nos,2*n_nos); % Vn*-w e Mn*-dwdn por no
G1 = zeros(2*n_nos,2*n_nos); % w*-Vn e dwdn*-Mn por no

for I = 1 : n_nos % Percorre as Fontes colocadas nos Elementos
    
    % Coordenadas do ponto fonte
    X0 = NOS(I,2:3);
    
    % Descobre qual elemento contem o ponto fonte
    for J = 1 : size(ELGEO,1)
        if I>=3*J-2 && I<=3*J
            elemfonte=J;
        end
    end
    noi=3*elemfonte-2; nof=3*elemfonte; %Nos inis e fins do elemfonte
    xi1=NOS(noi,2); yi1=NOS(noi,3); xi3=NOS(nof,2); yi3=NOS(nof,3);
    % Calculo da normal do ponto fonte a partir do elemento
    Li=sqrt((xi3-xi1)^2+(yi3-yi1)^2); % L do Elemfonte
    seno=(yi3-yi1)/Li; coss=(xi3-xi1)/Li; % Componentes do vetor tangente
    Ni=[seno;-coss]; % Vetor normal ao ponto fonte
    
    for J = 1 : n_el % Percorre os elementos, H1G1
        % Fontes nos nohs, integracao nos elementos
        
        % Coordenadas dos geo do elemento
        geo1 = ELGEO(J,2); X1 = GEO(geo1,2:3);
        geo2 = ELGEO(J,3); X2 = GEO(geo2,2:3);
        geo3 = ELGEO(J,4); X3 = GEO(geo3,2:3);
                
        % Calculo de h1 e g1 para cada elemento
        % Integracao nao-singular (regular)
        [h1,g1] = n3_h1g1Regular(X0,X1,X2,X3,Ni,D,v); %(2x6)
        
        % Integracao singular, teste fonte==dominio
        if I==3*J-2 % A fonte eh o noh 1 do elem J
            ksi = -2/3; col1 = 1; col2 = 2;
            
        elseif I==3*J-1 % A fonte eh o noh 2 do elem J
            ksi = 0; col1 = 3; col2 = 4;
                    
        elseif I==3*J % A fonte eh o noh 3 do elem J
            ksi = 2/3; col1 = 5; col2 = 6;
        end
        
        if I==3*J-2 || I==3*J-1 || I==3*J % If singular
            [hs,gs,HH]=n3_h1g1Sing5(X1,X2,X3,ksi,D,v); %(2x2,2x2,1x3)
            h1(:,col1:col2)=hs; % associa os 4 valores (2x2)
            g1(1,col2)=gs(1,2); % nao usa gs11
            g1(2,col1)=gs(2,1); % nao usa gs11
            g1(2,col2)=gs(2,2); % nao usa gs11
            h1(2,1)=HH(1); % associa 3 Hadamard (1x3)
            h1(2,3)=HH(2); % associa 3 Hadamard (1x3)
            h1(2,5)=HH(3); % associa 3 Hadamard (1x3)
            h1(1,2)=hMnA(2,X1,X3,v,ksi); %analit.
            h1(1,4)=hMnA(4,X1,X3,v,ksi); %analit.
            h1(1,6)=hMnA(6,X1,X3,v,ksi); %analit.
            g1(1,1)=gA(1,X1,X3,D,ksi); %analit.
            g1(2,2)=gA(2,X1,X3,D,ksi); %analit.
            g1(1,3)=gA(3,X1,X3,D,ksi); %analit.
            g1(2,4)=gA(4,X1,X3,D,ksi); %analit.
            g1(1,5)=gA(5,X1,X3,D,ksi); %analit.
            g1(2,6)=gA(6,X1,X3,D,ksi); %analit.
        end
        
        % Organiza componentes bloco por bloco (2x6)
        % junta as matrizes de elemento para formar matriz global
        H1(2*I-1:2*I,6*J-5:6*J)=h1; %contem Vn, dVn/dni, Mn e dMn/dni
        G1(2*I-1:2*I,6*J-5:6*J)=g1; %contem w, dw/dni, dw/dn e d2w/dndni
    end
end

H1=H1+1/2*eye(size(H1));