function [H2,G2]=n2montaH2G2sing3(CANTOS,ELGEO,GEO,D,v)

n_el = size(ELGEO,1); % Numero total de elementos
n_cantos = size(CANTOS,1); % Numero total de cantos

H2 = zeros(n_cantos,6*n_el);
G2 = zeros(n_cantos,6*n_el);

for I = 1 : n_cantos % Percorre as Fontes colocadas nos Cantos
    
    % Coordenadas do Ponto Fonte (Canto)
    geocanto=CANTOS(I,2);
    Xci=GEO(geocanto,2:3);
    elemAntes=CANTOS(I,3);
    elemDepois=CANTOS(I,4);
        
    for J = 1 : n_el % Percorre os elementos, H2G2
        % Fontes nos cantos, integracao nos elementos
                
        % Coordenadas dos tres geo do elemento
        geo1 = ELGEO(J,2); X1 = GEO(geo1,2:3);
        geo2 = ELGEO(J,3); X2 = GEO(geo2,2:3);
        geo3 = ELGEO(J,4); X3 = GEO(geo3,2:3);
                     
        % Integracao Regular
        [h2,g2] = n3_h2g2Regular(Xci,X1,X2,X3,D,v);
        
        % Integracao singular (apenas Mn, colunas pares)
        % teste elemento_de_integracao contem canto-fonte?
        if J==elemAntes || J==elemDepois
            if J==elemAntes % J eh Elem Antes do canto I
                h2cs=n3_h2CantoSing(X1,X2,X3,v,1); %colunas pares
            elseif J==elemDepois % J eh Elem Depois do canto I
                h2cs=n3_h2CantoSing(X1,X2,X3,v,2); %colunas pares
            end
            h2(1,2)=h2cs(1,2);h2(1,4)=h2cs(1,4);h2(1,6)=h2cs(1,6);
        end
               
        % Organiza componentes bloco por bloco (1x6)
        H2(I,6*J-5:6*J)=h2(:); % contem Vn e Mn
        G2(I,6*J-5:6*J)=g2(:); % contem w e dw/dn
    end
end