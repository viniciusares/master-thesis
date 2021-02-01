function C=CargasDominio(NOS,ELGEO,Pcels,Celulas,D,tipo)
% 2015/06/08 Integracao numerica bidimensional linear(4pontos)
% integra todas celulas para todas fontes
% tipo 2: w e dwdni (valido para nos do contorno) "2 por noh"
% tipo 1: w         (valido para cantos e p_ints) "1 por noh"

% 4 Pontos de Gauss Regulares
ksi=[-0.861136311594053 -0.339981043584856 ...
    0.339981043584856 0.861136311594053];
% 4 Pesos de Gauss Regulares
pg=[0.347854845137454 0.652145154862546 ...
    0.652145154862546 0.347854845137454];

npg=length(ksi); % Numero de pontos de Gauss
n_nos=size(NOS,1); %vale para nos, cantos, ou p_ints
C=zeros(tipo*n_nos,1); % tipo (1 por noh) ou (2 por noh)

for I=1:n_nos %percorre nos
    
    % Dados do ponto fonte
    x0=NOS(I,2); y0=NOS(I,3); % ponto fonte
    if tipo==2
        % Descobre qual elemento contem o ponto fonte
        for J = 1 : size(ELGEO,1)
            if I>=3*J-2 && I<=3*J
                elemfonte=J;
            end
        end
        noi=3*elemfonte-2; nof=3*elemfonte; %Nos inis e fins do elemfonte
        xi1=NOS(noi,2); yi1=NOS(noi,3); xi3=NOS(nof,2); yi3=NOS(nof,3);
        Li=sqrt((xi3-xi1)^2+(yi3-yi1)^2); % L do Elemfonte
        Ni=[(yi3-yi1)/Li;-(xi3-xi1)/Li]; % Vetor normal ao ponto fonte
    else
        Ni=[1;1]; %nao usado
    end
    
    p=zeros(1,4); % inicializacao
    if tipo==2          
        q=zeros(1,4);   
    end                 
    
    for J=1:size(Celulas,1) %percorre celulas
        InciCelulaJ=Celulas(J,2:5); %pontos que compoe a celula J
        ponto1=InciCelulaJ(1);
        ponto2=InciCelulaJ(2);
        ponto3=InciCelulaJ(3);
        ponto4=InciCelulaJ(4);
        xyc1=Pcels(ponto1,2:4); %coords xy e carga para ponto 1
        xyc2=Pcels(ponto2,2:4); %coords xy e carga para ponto 2
        xyc3=Pcels(ponto3,2:4); %coords xy e carga para ponto 3
        xyc4=Pcels(ponto4,2:4); %coords xy e carga para ponto 4
        xalfa=[xyc1(1);xyc2(1);xyc3(1);xyc4(1)]; % 4 valores de x
        yalfa=[xyc1(2);xyc2(2);xyc3(2);xyc4(2)]; % 4 valores de y
        cxy=[xyc1(3);xyc2(3);xyc3(3);xyc4(3)]; % 4 valores de cargas
      
        for K=1:npg % percorre eta1
            for M=1:npg % percorre eta2

                eta1=ksi(K); pg1=pg(K); %eta1 e eta2 sao pontos de Gauss
                eta2=ksi(M); pg2=pg(M); %pg1 e pg2 sao pesos de Gauss
                
                % Funcoes de forma usadas em x, y, p e q
                N1=1/4*(1+eta1)*(1+eta2); 
                N2=1/4*(1-eta1)*(1+eta2); 
                N3=1/4*(1-eta1)*(1-eta2);
                N4=1/4*(1+eta1)*(1-eta2);
                Nalfa=[N1 N2 N3 N4];
                
                % Para JAC
                dN1deta1=1/4*(1+eta2); dN1deta2=1/4*(1+eta1);
                dN2deta1=-1/4*(1+eta2); dN2deta2=1/4*(1-eta1);
                dN3deta1=-1/4*(1-eta2); dN3deta2=-1/4*(1-eta1);
                dN4deta1=1/4*(1-eta2); dN4deta2=-1/4*(1+eta1);
                dNalfadeta1=[dN1deta1 dN2deta1 dN3deta1 dN4deta1];
                dNalfadeta2=[dN1deta2 dN2deta2 dN3deta2 dN4deta2];
                
                dxdeta1=dNalfadeta1*xalfa; dxdeta2=dNalfadeta2*xalfa;
                dydeta1=dNalfadeta1*yalfa; dydeta2=dNalfadeta2*yalfa;
                JAC=det([dxdeta1 dydeta1;dxdeta2 dydeta2]);

                % Calcula as coordenadas do ponto de integracao 
                x=Nalfa*xalfa;
                y=Nalfa*yalfa;
                
                % Parametros               
                X=[x,y]; X0=[x0,y0]; % R diferente de 0, nao singular
                
                v=1; %nao usado
                [~,~,w,dwdn]=n4nucleos(X0,X,Ni,D,v);
                
                p=p+w*Nalfa*JAC*pg1*pg2;
                
                if tipo==2
                    dwdni=-dwdn;
                    q=q+dwdni*Nalfa*JAC*pg1*pg2;
                end                              
            end
        end
    end
    
    if tipo==1
        C(I,1)=p*cxy;
    elseif tipo==2
        C(2*I-1,1)=p*cxy;
        C(2*I,1)=q*cxy;
    else
        dips('ERRO:tipo so pode ser (1 por noh) ou (2 por noh)')
    end
    
end
