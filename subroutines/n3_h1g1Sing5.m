function [h,g,HH]=n3_h1g1Sing5(X1,X2,X3,ksi0,D,v)

x1=X1(1);y1=X1(2);
x2=X2(1);y2=X2(2);
x3=X3(1);y3=X3(2);

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
if(npgr~=npgs) disp('ERRO:nGauss_REG~=nGauss_SING')
else N=npgr;
end

% Singularidade Logaritmica
% fonte pertence ao elemento como noh ou como canto
% MnSing=-(1+v)/(4*pi)*(1+log(R))+(1-ni)/(8*pi);
%OK com kernels01marco15 30/04/15 %Difere do Regular
% MnSing contem parte singular+regular
MnReg=-(1+v)/(4*pi)+(1-v)/(8*pi); %parte regular do Mn singular

% d2wdnidnSing=-1/(8*pi*D)*(2*log(R)+1);
%OK com kernels01marco15 30/04/15 %Difere do Regular
% d2wdnidnSing contem parte singular+regular
d2wdndniReg=-1/(8*pi*D); %parte regular do d2wdnidn singular

% Para singularidade logaritmica
hs1=0; hs2=0;           % Inicializa a integral singular
hns1=0; hns2=0; hns3=0; % Inicializa a integral nao singular
gs1=0; gs2=0;           % Inicializa a integral singular
gns1=0; gns2=0; gns3=0; % Inicializa a integral nao singular
hreg=0; greg=0;         % Inicializa integral do termo regular

% Calculo das Singularidades Logaritmicas
if(ksi0==-2/3) % O ponto fonte eh o primeiro
    for I=1:N % Loop sobre os Pontos de Gauss (integracao numerica)
        
        % Integracao Regular (padrao) Mn e d2wdnidn
        RB=sqrt(((x1-2*x2+x3)*(ksi(I)+ksi0)+x3-x1)^2+...
                ((y1-2*y2+y3)*(ksi(I)+ksi0)+y3-y1)^2);
        qsi1=-1/6*(ksi(I)+5);
        qsi2=1/6*(5*ksi(I)+1);
        hns1=hns1-1*(-(1+v)/(4*pi))*(log(RB)-log(2))*Nd(ksi(I),1)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        gns1=gns1-(-1/(4*pi*D)*(log(RB)-log(2))*Nd(ksi(I),1)*...
            JACxy(X1,X2,X3,ksi(I),0))*pg(I);
        hns2=hns2-(1+v)/(24*pi)*log(3)*Nd(qsi1,1)*...
            JACxy(X1,X2,X3,qsi1,0)*pg(I);
        gns2=gns2-1/(24*pi*D)*log(3)*Nd(qsi1,1)*...
            JACxy(X1,X2,X3,qsi1,0)*pg(I);
        hns3=hns3-5*(1+v)/(24*pi)*log(3/5)*Nd(qsi2,1)*...
            JACxy(X1,X2,X3,qsi2,0)*pg(I);
        gns3=gns3-5/(24*pi*D)*log(3/5)*Nd(qsi2,1)*...
            JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        % Integracao Singular (logaritma) Mn e d2wdnidn
        qsi1l=-1/3*(eta(I)+2);
        qsi2l=1/3*(5*eta(I)-2);
        hs1=hs1-(1+v)/(12*pi)*Nd(qsi1l,1)*...
            JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        gs1=gs1-1/(12*pi*D)*Nd(qsi1l,1)*...
            JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        hs2=hs2-5*(1+v)/(12*pi)*Nd(qsi2l,1)*...
            JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        gs2=gs2-5/(12*pi*D)*Nd(qsi2l,1)*...
            JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
        % Parte regular
        hreg=hreg-MnReg*Nd(ksi(I),1)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        greg=greg-d2wdndniReg*Nd(ksi(I),1)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
    end
end
if(ksi0==0) % O ponto fonte eh o segundo
    for I=1:N
        
        % Integracao Regular (padrao)
        RB=sqrt(((x1-2*x2+x3)*(ksi(I)+ksi0)+x3-x1)^2+...
                ((y1-2*y2+y3)*(ksi(I)+ksi0)+y3-y1)^2);
        hns1=hns1-1*(-(1+v)/(4*pi))*((log(RB)-log(2))*Nd(ksi(I),2)*...
            JACxy(X1,X2,X3,ksi(I),0))*pg(I);
        gns1=gns1-(-1/(4*pi*D))*((log(RB)-log(2))*Nd(ksi(I),2)*...
            JACxy(X1,X2,X3,ksi(I),0))*pg(I);
        hns2=0; hns3=0; gns2=0; gns3=0;
        
        % Integracao Singular (logaritma)
        qsil=eta(I);
        hs1=hs1-(1+v)/(2*pi)*Nd(qsil,2)*...
            JACxy(X1,X2,X3,qsil,0)*rho(I);
        gs1=gs1-1/(2*pi*D)*Nd(qsil,2)*...
            JACxy(X1,X2,X3,qsil,0)*rho(I);
        hs2=0; gs2=0;
        
        % Parte regular
        hreg=hreg-MnReg*Nd(ksi(I),2)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        greg=greg-d2wdndniReg*Nd(ksi(I),2)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
    end
end
if(ksi0==2/3) % O ponto fonte eh o terceiro
    for I=1:N
        
        % Integracao Regular (padrao)
        RB=sqrt(((x1-2*x2+x3)*(ksi(I)+ksi0)+x3-x1)^2+...
                ((y1-2*y2+y3)*(ksi(I)+ksi0)+y3-y1)^2);
        qsi1=-1/6*(5*ksi(I)+1);
        qsi2=1/6*(ksi(I)+5);
        hns1=hns1-1*(-(1+v)/(4*pi))*(log(RB)-log(2))*Nd(ksi(I),3)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        gns1=gns1-(-1/(4*pi*D))*(log(RB)-log(2))*Nd(ksi(I),3)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        hns2=hns2-5*(1+v)/(24*pi)*log(3/5)*Nd(qsi1,3)*...
            JACxy(X1,X2,X3,qsi1,0)*pg(I);
        gns2=gns2-5/(24*pi*D)*log(3/5)*Nd(qsi1,3)*...
            JACxy(X1,X2,X3,qsi1,0)*pg(I);
        hns3=hns3-(1+v)/(24*pi)*log(3)*Nd(qsi2,3)*...
            JACxy(X1,X2,X3,qsi2,0)*pg(I);
        gns3=gns3-1/(24*pi*D)*log(3)*Nd(qsi2,3)*...
            JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        % Integracao Singular (logaritma)
        qsi1l=-1/3*(5*eta(I)-2);
        qsi2l=1/3*(eta(I)+2);
        hs1=hs1-5*(1+v)/(12*pi)*Nd(qsi1l,3)*...
            JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        gs1=gs1-5/(12*pi*D)*Nd(qsi1l,3)*...
            JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        hs2=hs2-(1+v)/(12*pi)*Nd(qsi2l,3)*...
            JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        gs2=gs2-1/(12*pi*D)*Nd(qsi2l,3)*...
            JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
        % Parte regular
        hreg=hreg-MnReg*Nd(ksi(I),3)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        greg=greg-d2wdndniReg*Nd(ksi(I),3)*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
    end
end

% Vn
h(1,1)=0; % Eh zero de fato

% Mn Singularidade Logaritmica
hns=hns1+hns2+hns3;
hs=hs1+hs2;
h(1,2)=hreg+hs+hns;

% dVndni HipersingularHadamard dVndni=(1+ni)/(4*pi*R^2);
h(2,1)=0; % pois sera sobrescrito por HH em n2montaHG

ep=2e-15; % tolerancia de comparacao
ruidoX=x2-(x1+x3)/2; ruidoY=y2-(y1+y3)/2;
if abs(ruidoX)>ep || abs(ruidoY)>ep
    disp('ERROn3_h1g1Sing5:Foi encontrado um elemento curvo')
    ruidoX
    ruidoY
else
    L=sqrt((x3-x1)^2+(y3-y1)^2);
end
ksil=ksi0;
HH(1)=3*(1+v)/(8*pi*L)*((3*ksil-1)*log(abs((1-ksil)/(1+ksil)))+...
    (6*ksil^2-2*ksil-3)/(ksil^2-1));
HH(2)=(1+v)/(4*pi*L)*(9*ksil*log(abs((1+ksil)/(1-ksil)))-...
    (18*ksil^2-13)/(ksil^2-1));
HH(3)=3*(1+v)/(8*pi*L)*((3*ksil+1)*log(abs((1-ksil)/(1+ksil)))+...
    (6*ksil^2+2*ksil-3)/(ksil^2-1));

% dMndni
h(2,2)=0; % Eh zero de fato

% w NaoSingular
g(1,1)=0; %valor nao atribuido em n2montaHG
%portanto nao sobrescreve n3calc_hgRegular

% dwdn
g(1,2)=0; % Eh zero de fato

% dwdni
g(2,1)=0; % Eh zero de fato

% d2wdndni Singularidade Logaritmica
gns=gns1+gns2+gns3;
gs=gs1+gs2;
g(2,2)=greg+gs+gns;

