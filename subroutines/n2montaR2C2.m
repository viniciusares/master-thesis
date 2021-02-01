function [R2,C2]=n2montaR2C2(CANTOS,ELGEO,GEO,D,v,H2)

n_cantos = size(CANTOS,1); % Numero total de cantos

R2 = zeros(n_cantos,n_cantos);
C2 = zeros(n_cantos,n_cantos);

for I = 1 : n_cantos % Percorre as Fontes colocadas nos Cantos
    
    % Coordenadas do Ponto Fonte (Canto)
    geocanto=CANTOS(I,2);
    Xci=GEO(geocanto,2:3);
    
    for J = 1 : n_cantos % Percorre os cantos, R2C2
        % Fontes nos cantos, avaliacao nos cantos
        
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
        
        % Avaliacao Regular (I~=J)
        [r2,c2] = n3_r2c2Regular(Xci,Xcj,Nj,NJ,Sj,SJ,D,v);
        
        % r2 e c2 Singulares
        if I==J
            r2=0; % corpo_rigido, sera sobrescrito abaixo
            c2=0; % R^2 mais forte que log(R), tende a zero
        end
                
        % Organiza componentes escalar por escalar (1x1)
        R2(I,J)=r2; % contem Rc apenas
        C2(I,J)=c2; % contem wc apenas
    end
    
    R2(I,I)=-sum(R2(I,:))-sum(H2(I,1:2:end)); % Corpo_Rigido    
end
