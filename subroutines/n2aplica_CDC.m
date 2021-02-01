function [A,b]=n2aplica_CDC(CDC,Hloc,Gloc,CTOS)

u=zeros(1,6*size(CDC,1)+size(CTOS,1));

for I=1:size(CDC,1) % Percorre os elementos
    
    if CDC(I,2)==1 % engastado, troca 2 (ambas)
        aux=-Hloc(:,6*I-5:6*I);
        Hloc(:,6*I-5:6*I)=-Gloc(:,6*I-5:6*I);
        Gloc(:,6*I-5:6*I)=aux;
                                          
    elseif CDC(I,2)==2 % apoiado,   troca 1 (1as colunas)
        aux=-Hloc(:,6*I-5:2:6*I-1); % armazena 3 col de H que serao trocadas
        Hloc(:,6*I-5:2:6*I-1)=-Gloc(:,6*I-5:2:6*I-1); % H recebe 3 cols de G
        Gloc(:,6*I-5:2:6*I-1)=aux; % G recebe 3 cols de H que estavam em aux
              
    elseif CDC(I,2)==3 % guia,      troca 1 (2as colunas)
        aux=-Hloc(:,6*I-4:2:6*I); % armazena 3 col de H que serao trocadas
        Hloc(:,6*I-4:2:6*I)=-Gloc(:,6*I-4:2:6*I); % H recebe 3 cols de G
        Gloc(:,6*I-4:2:6*I)=aux; % G recebe 3 cols de H que estavam em aux
                        
    elseif CDC(I,2)==4 % livre,     nao troca nada
        
    
    else disp('ERRO:tipoCDC deve estar entre 1 e 4')
    end
    
    u(6*I-5:6*I)=CDC(I,3:8); % armazena valores CDC em vetor coluna
end

for I=1:size(CTOS,1) % percorre os cantos
    J=I+6*size(CDC,1); % colunas de cantos, apos elementos
    if CTOS(I,5) == 1 %canto preso-1:troca, else=solto-2:nao faz nada
        aux=-Hloc(:,J);         % armazena 1 col para troca
        Hloc(:,J)=-Gloc(:,J);   % troca 1 col
        Gloc(:,J)=aux;          % troca 1 col (recuperando)
    end
    
    u(J)=CTOS(I,6); % armazena valores CDCcantos no final do vetor
end

b=Gloc*u'; % parte correspondente ao que se conhece
A=Hloc; % matriz correspondente aas incognitas
