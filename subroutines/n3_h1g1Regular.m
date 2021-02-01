function [h,g]=n3_h1g1Regular(X0,X1,X2,X3,Ni,D,v)

x1=X1(1);y1=X1(2);
x2=X2(1);y2=X2(2);
x3=X3(1);y3=X3(2);

% Pontos de Gauss
ksi = [-0.973906528517172 0.973906528517172 -0.865063366688985 0.865063366688985 ...
        -0.679409568299024 0.679409568299024 -0.433395394129247 0.433395394129247 ...
        -0.148874338981631 0.148874338981631];
% Pesos de Gauss
pg = [ 0.066671344308688 0.066671344308688 0.149451349150581 0.149451349150581 ...
        0.219086362515982 0.219086362515982 0.269266719309996 0.269266719309996 ...
        0.295524224714753 0.295524224714753];
    
npg=length(ksi); % Numero de pontos de Gauss

h=zeros(2,6); % Inicializa a matriz h do elemento
g=zeros(2,6); % Inicializa a matriz g do elemento

for I=1:npg % Percorre os pontos de Gauss
        
    % Calcula as coordenadas do ponto de integracao
    x=Nc(ksi(I),1)*x1+Nc(ksi(I),2)*x2+Nc(ksi(I),3)*x3;
    y=Nc(ksi(I),1)*y1+Nc(ksi(I),2)*y2+Nc(ksi(I),3)*y3; 
    X=[x,y];
    
    % Componentes x e y do vetor tangente ao elemento
    sx=1/JACxy(X1,X2,X3,ksi(I),0)*JACxy(X1,X2,X3,ksi(I),1); 
    sy=1/JACxy(X1,X2,X3,ksi(I),0)*JACxy(X1,X2,X3,ksi(I),2);
    nx=sy; ny=-sx; N=[nx;ny]; % Vetor normal ao elemento
    
    % Solucoes Fundamentais (nucleos)
    [Vn,Mn,w,dwdn]=n4nucleos(X0,X,N,D,v);
    [dVndni,dMndni,dwdni,d2wdndni]=n4nucleosH(X0,X,Ni,N,D,v);
    
    JAC=JACxy(X1,X2,X3,ksi(I),0);
    %-------Integracao componentes h------------------------
    h(1,1)=h(1,1)+Vn*Nd(ksi(I),1)*JAC*pg(I);
    h(1,3)=h(1,3)+Vn*Nd(ksi(I),2)*JAC*pg(I);
    h(1,5)=h(1,5)+Vn*Nd(ksi(I),3)*JAC*pg(I);
      
    h(2,1)=h(2,1)+dVndni*Nd(ksi(I),1)*JAC*pg(I);
    h(2,3)=h(2,3)+dVndni*Nd(ksi(I),2)*JAC*pg(I);
    h(2,5)=h(2,5)+dVndni*Nd(ksi(I),3)*JAC*pg(I);
    
    h(1,2)=h(1,2)-Mn*Nd(ksi(I),1)*JAC*pg(I);
    h(1,4)=h(1,4)-Mn*Nd(ksi(I),2)*JAC*pg(I);
    h(1,6)=h(1,6)-Mn*Nd(ksi(I),3)*JAC*pg(I);
    
    h(2,2)=h(2,2)-dMndni*Nd(ksi(I),1)*JAC*pg(I);
    h(2,4)=h(2,4)-dMndni*Nd(ksi(I),2)*JAC*pg(I);
    h(2,6)=h(2,6)-dMndni*Nd(ksi(I),3)*JAC*pg(I);

    %-------Integracao componentes g------------------------
    g(1,1)=g(1,1)+w*Nd(ksi(I),1)*JAC*pg(I);
    g(1,3)=g(1,3)+w*Nd(ksi(I),2)*JAC*pg(I);
    g(1,5)=g(1,5)+w*Nd(ksi(I),3)*JAC*pg(I);
       
    g(2,1)=g(2,1)+dwdni*Nd(ksi(I),1)*JAC*pg(I);
    g(2,3)=g(2,3)+dwdni*Nd(ksi(I),2)*JAC*pg(I);
    g(2,5)=g(2,5)+dwdni*Nd(ksi(I),3)*JAC*pg(I);
    
    g(1,2)=g(1,2)-dwdn*Nd(ksi(I),1)*JAC*pg(I);
    g(1,4)=g(1,4)-dwdn*Nd(ksi(I),2)*JAC*pg(I);
    g(1,6)=g(1,6)-dwdn*Nd(ksi(I),3)*JAC*pg(I);
    
    g(2,2)=g(2,2)-d2wdndni*Nd(ksi(I),1)*JAC*pg(I);
    g(2,4)=g(2,4)-d2wdndni*Nd(ksi(I),2)*JAC*pg(I);
    g(2,6)=g(2,6)-d2wdndni*Nd(ksi(I),3)*JAC*pg(I);

end