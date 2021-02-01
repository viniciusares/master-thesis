function [h2,g2]=n3_h2g2Regular(X0,X1,X2,X3,D,v)

x1=X1(1);y1=X1(2);
x2=X2(1);y2=X2(2);
x3=X3(1);y3=X3(2);

% Pontos de Gauss
ksi = [-0.97390652851717 0.97390652851717 -0.86506336668898 0.86506336668898 ...
        -0.67940956829902 0.67940956829902 -0.43339539412925 0.43339539412925 ...
        -0.14887433898163 0.14887433898163];
% Pesos de Gauss
pg = [ 0.06667134430868 0.06667134430868 0.14945134915058 0.14945134915058 ...
        0.21908636251598 0.21908636251598 0.26926671931000 0.26926671931000 ...
        0.29552422471475 0.29552422471475];
    
npg=length(ksi); % Numero de pontos de Gauss

h2=zeros(1,6); % Inicializa a matriz h2 do elemento
g2=zeros(1,6); % Inicializa a matriz g2 do elemento

for I=1:npg % Percorre os pontos de Gauss

    % Calcula as coordenadas do ponto de integracao
    x=Nc(ksi(I),1)*x1+Nc(ksi(I),2)*x2+Nc(ksi(I),3)*x3;
    y=Nc(ksi(I),1)*y1+Nc(ksi(I),2)*y2+Nc(ksi(I),3)*y3; 
    X=[x,y];
    
    % Componentes x e y do vetor tangente (do ELEMENTO)
    sx=1/JACxy(X1,X2,X3,ksi(I),0)*JACxy(X1,X2,X3,ksi(I),1); 
    sy=1/JACxy(X1,X2,X3,ksi(I),0)*JACxy(X1,X2,X3,ksi(I),2);
    nx=sy; ny=-sx; N=[nx;ny]; % Vetor normal

    [Vn,Mn,w,dwdn]=n4nucleos(X0,X,N,D,v);

    JAC=JACxy(X1,X2,X3,ksi(I),0);
    %-------Integracao componentes h------------------------
    h2(1,1)=h2(1,1)+Vn*Nd(ksi(I),1)*JAC*pg(I);
    h2(1,3)=h2(1,3)+Vn*Nd(ksi(I),2)*JAC*pg(I);
    h2(1,5)=h2(1,5)+Vn*Nd(ksi(I),3)*JAC*pg(I);
          
    h2(1,2)=h2(1,2)-Mn*Nd(ksi(I),1)*JAC*pg(I);
    h2(1,4)=h2(1,4)-Mn*Nd(ksi(I),2)*JAC*pg(I);
    h2(1,6)=h2(1,6)-Mn*Nd(ksi(I),3)*JAC*pg(I);
    
    %-------Integracao componentes g------------------------
    g2(1,1)=g2(1,1)+w*Nd(ksi(I),1)*JAC*pg(I);
    g2(1,3)=g2(1,3)+w*Nd(ksi(I),2)*JAC*pg(I);
    g2(1,5)=g2(1,5)+w*Nd(ksi(I),3)*JAC*pg(I);
           
    g2(1,2)=g2(1,2)-dwdn*Nd(ksi(I),1)*JAC*pg(I);
    g2(1,4)=g2(1,4)-dwdn*Nd(ksi(I),2)*JAC*pg(I);
    g2(1,6)=g2(1,6)-dwdn*Nd(ksi(I),3)*JAC*pg(I);

end