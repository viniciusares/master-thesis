function h2cs=n3_h2CantoSing(X1,X2,X3,v,AD)

% Pontos de Gauss Padrao
ksi = [-0.97390652851717 0.97390652851717 -0.86506336668898 0.86506336668898 ...
    -0.67940956829902 0.67940956829902 -0.43339539412925 0.43339539412925 ...
    -0.14887433898163 0.14887433898163];
% Pesos de Gauss Padrao
pg = [ 0.06667134430868 0.06667134430868 0.14945134915058 0.14945134915058 ...
    0.21908636251598 0.21908636251598 0.26926671931000 0.26926671931000 ...
    0.29552422471475 0.29552422471475];
% Pontos de Gauss logaritmico
eta = [0.90425944e-2 0.53971054e-1	0.13531134 0.24705169 0.38021171 ...
    0.52379159 0.66577472 0.79419019 0.89816102 0.96884798];
% Pesos de Gauss Logaritmico
rho = [0.12095474 0.18636310 0.19566066 0.17357723 0.13569597 0.93647084e-1 ...
    0.55787938e-1 0.27159893e-1 0.95151992e-2 0.16381586e-2];

npgr=length(ksi); % Numero de Pontos de Gauss Regular (padrao)
npgs=length(eta); % Numero de Pontos de Gauss Singular (logaritmico)
if(npgr~=npgs) disp('ERROn3_h2CantoSing:nGauss_REG~=nGauss_SING')
else N=npgr;
end

x1=X1(1); x2=X2(1); x3=X3(1); y1=X1(2); y2=X2(2); y3=X3(2);
ep=2e-15; % tolerancia de comparacao
ruidoX=x2-(x1+x3)/2; ruidoY=y2-(y1+y3)/2;
if abs(ruidoX)>ep || abs(ruidoY)>ep
    disp('ERROn3_h2CantoSing:Foi encontrado um elemento curvo')
    ruidoX
    ruidoY
else
    L=sqrt((x3-x1)^2+(y3-y1)^2);
end

h2reg=zeros(1,6);
h2sing=zeros(1,6);
% MnSing=-(1+v)/(4*pi)*log((1-ksi)/2);         % nao usado
% MnReg=-(1+v)/(4*pi)*(1+log(L))+(1-v)/(8*pi); % nao usado

for I=1:N % percorre os pontos de gauss
    
    % Integracao Regular (nao muda se antes ou depois)
    h2reg(1,2)=h2reg(1,2)+(L*(1+v)/(8*pi)*(1+log(L))...
        -L*(1-v)/(16*pi))*Nd(ksi(I),1)*pg(I);
    
    h2reg(1,4)=h2reg(1,4)+(L*(1+v)/(8*pi)*(1+log(L))...
        -L*(1-v)/(16*pi))*Nd(ksi(I),2)*pg(I);
    
    h2reg(1,6)=h2reg(1,6)+(L*(1+v)/(8*pi)*(1+log(L))...
        -L*(1-v)/(16*pi))*Nd(ksi(I),3)*pg(I);
        
    % Definicao Antes-Depois para Integracao Singular
    switch AD
        case 1 % Antes
            %ETAil=(1-eta(I))/2;
            %fnil1=Nd(1-2*ETAil,1);
            %fnil2=Nd(1-2*ETAil,2);
            %fnil3=Nd(1-2*ETAil,3);
            fnil1=Nd(1-2*eta(I),1);
            fnil2=Nd(1-2*eta(I),2);
            fnil3=Nd(1-2*eta(I),3);
        case 2 % Depois
            %ETAil=(1+eta(I))/2;
            %fnil1=Nd(2*ETAil-1,1);
            %fnil2=Nd(2*ETAil-1,2);
            %fnil3=Nd(2*ETAil-1,3);
            fnil1=Nd(2*eta(I)-1,1);
            fnil2=Nd(2*eta(I)-1,2);
            fnil3=Nd(2*eta(I)-1,3);
    end
    % Integracao Singular    
    h2sing(1,2)=h2sing(1,2)-L*(1+v)/(4*pi)*fnil1*rho(I);    
    h2sing(1,4)=h2sing(1,4)-L*(1+v)/(4*pi)*fnil2*rho(I);   
    h2sing(1,6)=h2sing(1,6)-L*(1+v)/(4*pi)*fnil3*rho(I);
    
    %--------------------------------------------------
    h2cs=h2sing+h2reg;
    
end
