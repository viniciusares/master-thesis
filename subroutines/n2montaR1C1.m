function [R1,C1]=n2montaR1C1(NOS,ELGEO,GEO,CANTOS,D,v)

n_nos = size(NOS,1); % Numero total de nos
n_cantos = size(CANTOS,1); % Numero total de cantos

R1 = zeros(2*n_nos,n_cantos);
C1 = zeros(2*n_nos,n_cantos);

for I = 1 : n_nos % Percorre as Fontes colocadas nos Elementos
    
    % Coordenadas do ponto fonte
    X0 = NOS(I,2:3);
    
    % Descobre qual elemento corresponde ao ponto fonte
    for J = 1 : size(ELGEO,1)
        if I>=3*J-2 && I<=3*J
            elemfonte=J;
        end
    end
    noi=3*elemfonte-2; nof=3*elemfonte; %Nos inis e fins do elemfonte
    xi1=NOS(noi,2); yi1=NOS(noi,3); xi3=NOS(nof,2); yi3=NOS(nof,3);
    
    % Calculo da normal do ponto fonte a partir do elemento
    L=sqrt((xi3-xi1)^2+(yi3-yi1)^2); % L do Elemfonte
    seno=(yi3-yi1)/L; coss=(xi3-xi1)/L; % Componentes do vetor tangente
    Ni=[seno;-coss]; % Vetor normal ao ponto fonte
    
    for J = 1 : n_cantos % Percorre os cantos, R1C1
        % Fontes nos nohs, avaliacao nos cantos
        
        % Descobre coordenadas do canto
        geocanto=CANTOS(J,2);
        Xcj=GEO(geocanto,2:3);
        
        % Descobre coordenadas dos elems antes e depois do canto
        elemantes=CANTOS(J,3);
        elemdepois=CANTOS(J,4);
        
        GEOSantes(1:3)=ELGEO(elemantes,2:4);
        GEOSdepois(1:3)=ELGEO(elemdepois,2:4);
        
        xj3=GEO(GEOSantes(3),2); yj3=GEO(GEOSantes(3),3);
        xj1=GEO(GEOSantes(1),2); yj1=GEO(GEOSantes(1),3);
        xJ3=GEO(GEOSdepois(3),2); yJ3=GEO(GEOSdepois(3),3);
        xJ1=GEO(GEOSdepois(1),2); yJ1=GEO(GEOSdepois(1),3);
        
        % Calculo dos vetores antes e depois do canto
        La=sqrt((xj3-xj1)^2+(yj3-yj1)^2); % L do elem antes
        Ld=sqrt((xJ3-xJ1)^2+(yJ3-yJ1)^2); % L do elem depois
        seno=(yj3-yj1)/La; coss=(xj3-xj1)/La; % antes
        SENO=(yJ3-yJ1)/Ld; COSS=(xJ3-xJ1)/Ld; % depois
        % Vetores coluna para produto escalar
        Sj=[coss;seno]; Nj=[seno;-coss]; % tangente e normal antes
        SJ=[COSS;SENO]; NJ=[SENO;-COSS]; % tangente e normal depois
        
        % Calculo das submatrizes r1_el e c1_el (apenas nao sing)
        % Avaliacao Regular
        [r1,c1] = n3_r1c1Regular(X0,Xcj,Ni,Nj,NJ,Sj,SJ,D,v);
        % Nao ha Avaliacao Singular
        
        % Organiza componentes vetor por vetor (2x1)
        R1(2*I-1:2*I,J)=r1; %contem Rc e dRc/dni
        C1(2*I-1:2*I,J)=c1; %contem wc e dwc/dni
    end % fecha loop nos cantos-avaliacao
end % fecha loop nos nohs-fonte
