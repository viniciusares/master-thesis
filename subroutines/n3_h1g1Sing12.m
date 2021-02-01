function [hMn,HH,g11,g22]=n3_h1g1Sing12(X1,X2,X3,ksi0,D,v)

x1=X1(1);y1=X1(2);
x2=X2(1);y2=X2(2);
x3=X3(1);y3=X3(2);

[ksi,pg,eta,rho]=Gauss(10);

npgr=length(ksi); % Numero de Pontos de Gauss Regular (padrao)
npgs=length(eta); % Numero de Pontos de Gauss Singular (logaritmico)
if(npgr~=npgs) disp('ERROn2_h1g1Sing12:nGauss_REG~=nGauss_SING')
else N=npgr;
end

% Singularidade Logaritmica
% MnSing=-(1+v)/(4*pi)*(1+log(R))+(1-ni)/(8*pi);
% MnSing contem parte singular+regular
MnReg=-(1+v)/(4*pi)+(1-v)/(8*pi); %parte regular do Mn singular
% d2wdnidnSing=-1/(8*pi*D)*(2*log(R)+1);
% d2wdnidnSing contem parte singular+regular
d2wdndniReg=-1/(8*pi*D); %parte regular do d2wdnidn singular

% Para singularidade logaritmica, Inicializacao de somatorios
ini=[0 0 0];
hns1=ini; hns2=ini; hns3=ini; % Mn
hs1=ini; hs2=ini; hreg=ini;

g11ns1=ini; g11ns2=ini; g11ns3=ini; % w
g11s1=ini; g11s2=ini; g11reg=ini;

g22ns1=ini; g22ns2=ini; g22ns3=ini; % d2wdndi
g22s1=ini; g22s2=ini; g22reg=ini;

% Calculo das Singularidades Logaritmicas
for I=1:N % Loop sobre os Pontos de Gauss (integracao numerica)
    RB=sqrt(((x1-2*x2+x3)*(ksi(I)+ksi0)+x3-x1)^2+...
        ((y1-2*y2+y3)*(ksi(I)+ksi0)+y3-y1)^2);
    
    ihns1=-1*(-(1+v)/(4*pi))*(log(RB)-log(2))*...
        JACxy(X1,X2,X3,ksi(I),0)*pg(I);
    
    x0=x1*Nc(ksi0,1)+x2*Nc(ksi0,2)+x3*Nc(ksi0,3);
    y0=y1*Nc(ksi0,1)+y2*Nc(ksi0,2)+y3*Nc(ksi0,3);
    x=x1*Nc(ksi(I),1)+x2*Nc(ksi(I),2)+x3*Nc(ksi(I),3);
    y=y1*Nc(ksi(I),1)+y2*Nc(ksi(I),2)+y3*Nc(ksi(I),3);
    R=sqrt((x-x0)^2+(y-y0)^2);
    ig11ns1=-(1/(8*pi*D)*(R^2*log(RB)-R^2*log(2))*...
        JACxy(X1,X2,X3,ksi(I),0))*pg(I);
    
    ig22ns1=-(-1/(4*pi*D)*(log(RB)-log(2))*...
        JACxy(X1,X2,X3,ksi(I),0))*pg(I);
    
    ihreg=-MnReg*...
            JACxy(X1,X2,X3,ksi(I),0)*pg(I);
    ig22reg=-d2wdndniReg*...
        JACxy(X1,X2,X3,ksi(I),0)*pg(I);
        
    if(ksi0==-2/3) % O ponto fonte eh o primeiro
        qsi1=-1/6*(ksi(I)+5);
        qsi10=-1/6*(ksi0+5);
        qsi2=1/6*(5*ksi(I)+1);
        qsi20=1/6*(5*ksi0+1);
        qsi1l=-1/3*(eta(I)+2);
        eta0=0;                 % VERIFICAR!!!
        qsi1l0=-1/3*(eta0+2);
        qsi2l=1/3*(5*eta(I)-2);
        qsi2l0=1/3*(5*eta0-2);
        
        ihns2=-(1+v)/(24*pi)*log(3)*JACxy(X1,X2,X3,qsi1,0)*pg(I);
        
        x0=x1*Nc(qsi10,1)+x2*Nc(qsi10,2)+x3*Nc(qsi10,3);
        y0=y1*Nc(qsi10,1)+y2*Nc(qsi10,2)+y3*Nc(qsi10,3);
        x=x1*Nc(qsi1,1)+x2*Nc(qsi1,2)+x3*Nc(qsi1,3);
        y=y1*Nc(qsi1,1)+y2*Nc(qsi1,2)+y3*Nc(qsi1,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11ns2=1/(48*pi*D)*R^2*log(3)*JACxy(X1,X2,X3,qsi1,0)*pg(I);
        
        ig22ns2=-1/(24*pi*D)*log(3)*JACxy(X1,X2,X3,qsi1,0)*pg(I);
        ihns3=-5*(1+v)/(24*pi)*log(3/5)*JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        x0=x1*Nc(qsi20,1)+x2*Nc(qsi20,2)+x3*Nc(qsi20,3);
        y0=y1*Nc(qsi20,1)+y2*Nc(qsi20,2)+y3*Nc(qsi20,3);
        x=x1*Nc(qsi2,1)+x2*Nc(qsi2,2)+x3*Nc(qsi2,3);
        y=y1*Nc(qsi2,1)+y2*Nc(qsi2,2)+y3*Nc(qsi2,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11ns3=5/(48*pi*D)*R^2*log(3/5)*JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        ig22ns3=-5/(24*pi*D)*log(3/5)*JACxy(X1,X2,X3,qsi2,0)*pg(I);
        ihs1=-(1+v)/(12*pi)*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        x0=x1*Nc(qsi1l0,1)+x2*Nc(qsi1l0,2)+x3*Nc(qsi1l0,3);
        y0=y1*Nc(qsi1l0,1)+y2*Nc(qsi1l0,2)+y3*Nc(qsi1l0,3);
        x=x1*Nc(qsi1l,1)+x2*Nc(qsi1l,2)+x3*Nc(qsi1l,3);
        y=y1*Nc(qsi1l,1)+y2*Nc(qsi1l,2)+y3*Nc(qsi1l,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11s1=1/(24*pi*D)*R^2*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        ig22s1=-1/(12*pi*D)*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        ihs2=-5*(1+v)/(12*pi)*JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
        x0=x1*Nc(qsi2l0,1)+x2*Nc(qsi2l0,2)+x3*Nc(qsi2l0,3);
        y0=y1*Nc(qsi2l0,1)+y2*Nc(qsi2l0,2)+y3*Nc(qsi2l0,3);
        x=x1*Nc(qsi2l,1)+x2*Nc(qsi2l,2)+x3*Nc(qsi2l,3);
        y=y1*Nc(qsi2l,1)+y2*Nc(qsi2l,2)+y3*Nc(qsi2l,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11s2=5/(24*pi*D)*R^2*JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
        ig22s2=-5/(12*pi*D)*JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
    elseif(ksi0==0) % O ponto fonte eh o segundo
        qsi1=0;
        qsi2=0;
        qsi1l=eta(I);
        eta0=0;     % VERIFICAR!!!
        qsi1l0=eta0;
        qsi2l=0;
        ihns2=0;
        ig11ns2=0;
        ig22ns2=0;
        ihns3=0;
        ig11ns3=0;
        ig22ns3=0;
        ihs1=-(1+v)/(2*pi)*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        x0=x1*Nc(qsi1l0,1)+x2*Nc(qsi1l0,2)+x3*Nc(qsi1l0,3);
        y0=y1*Nc(qsi1l0,1)+y2*Nc(qsi1l0,2)+y3*Nc(qsi1l0,3);
        x=x1*Nc(qsi1l,1)+x2*Nc(qsi1l,2)+x3*Nc(qsi1l,3);
        y=y1*Nc(qsi1l,1)+y2*Nc(qsi1l,2)+y3*Nc(qsi1l,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11s1=1/(4*pi*D)*R^2*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        ig22s1=-1/(2*pi*D)*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        ihs2=0;
        ig11s2=0;
        ig22s2=0;
        
    elseif(ksi0==2/3) % O ponto fonte eh o terceiro
        qsi1=-1/6*(5*ksi(I)+1);
        qsi10=-1/6*(5*ksi0+1);
        qsi2=1/6*(ksi(I)+5);
        qsi20=1/6*(ksi0+5);
        qsi1l=-1/3*(5*eta(I)-2);
        eta0=0;                     % VERIFICAR !!!
        qsi1l0=-1/3*(5*eta0-2);
        qsi2l=1/3*(eta(I)+2);
        qsi2l0=1/3*(eta0+2);
        
        ihns2=-5*(1+v)/(24*pi)*log(3/5)*JACxy(X1,X2,X3,qsi1,0)*pg(I);
        
        x0=x1*Nc(qsi10,1)+x2*Nc(qsi10,2)+x3*Nc(qsi10,3);
        y0=y1*Nc(qsi10,1)+y2*Nc(qsi10,2)+y3*Nc(qsi10,3);
        x=x1*Nc(qsi1,1)+x2*Nc(qsi1,2)+x3*Nc(qsi1,3);
        y=y1*Nc(qsi1,1)+y2*Nc(qsi1,2)+y3*Nc(qsi1,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11ns2=5/(48*pi*D)*R^2*log(3/5)*JACxy(X1,X2,X3,qsi1,0)*pg(I);
        
        ig22ns2=-5/(24*pi*D)*log(3/5)*JACxy(X1,X2,X3,qsi1,0)*pg(I);
        
        ihns3=-(1+v)/(24*pi)*log(3)*JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        x0=x1*Nc(qsi20,1)+x2*Nc(qsi20,2)+x3*Nc(qsi20,3);
        y0=y1*Nc(qsi20,1)+y2*Nc(qsi20,2)+y3*Nc(qsi20,3);
        x=x1*Nc(qsi2,1)+x2*Nc(qsi2,2)+x3*Nc(qsi2,3);
        y=y1*Nc(qsi2,1)+y2*Nc(qsi2,2)+y3*Nc(qsi2,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11ns3=1/(48*pi*D)*R^2*log(3)*JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        ig22ns3=-1/(24*pi*D)*log(3)*JACxy(X1,X2,X3,qsi2,0)*pg(I);
        
        ihs1=-5*(1+v)/(12*pi)*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        x0=x1*Nc(qsi1l0,1)+x2*Nc(qsi1l0,2)+x3*Nc(qsi1l0,3);
        y0=y1*Nc(qsi1l0,1)+y2*Nc(qsi1l0,2)+y3*Nc(qsi1l0,3);
        x=x1*Nc(qsi1l,1)+x2*Nc(qsi1l,2)+x3*Nc(qsi1l,3);
        y=y1*Nc(qsi1l,1)+y2*Nc(qsi1l,2)+y3*Nc(qsi1l,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11s1=5/(24*pi*D)*R^2*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        ig22s1=-5/(12*pi*D)*JACxy(X1,X2,X3,qsi1l,0)*rho(I);
        
        ihs2=-(1+v)/(12*pi)*JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
        x0=x1*Nc(qsi2l0,1)+x2*Nc(qsi2l0,2)+x3*Nc(qsi2l0,3);
        y0=y1*Nc(qsi2l0,1)+y2*Nc(qsi2l0,2)+y3*Nc(qsi2l0,3);
        x=x1*Nc(qsi2l,1)+x2*Nc(qsi2l,2)+x3*Nc(qsi2l,3);
        y=y1*Nc(qsi2l,1)+y2*Nc(qsi2l,2)+y3*Nc(qsi2l,3);
        R=sqrt((x-x0)^2+(y-y0)^2);
        ig11s2=1/(24*pi*D)*R^2*JACxy(X1,X2,X3,qsi2l,0)*rho(I);
        
        ig22s2=-1/(12*pi*D)*JACxy(X1,X2,X3,qsi2l,0)*rho(I);
    end
        
    % Parte singular-regular: Mn, w e d2wdndni
        % ns1 com ksi(I)
    hns1(1)=hns1(1)+ihns1*Nd(ksi(I),1);
    hns1(2)=hns1(2)+ihns1*Nd(ksi(I),2);
    hns1(3)=hns1(3)+ihns1*Nd(ksi(I),3);

    g11ns1(1)=g11ns1(1)+ig11ns1*Nd(ksi(I),1);
    g11ns1(2)=g11ns1(2)+ig11ns1*Nd(ksi(I),2);
    g11ns1(3)=g11ns1(3)+ig11ns1*Nd(ksi(I),3);
    g22ns1(1)=g22ns1(1)+ig22ns1*Nd(ksi(I),1);
    g22ns1(2)=g22ns1(2)+ig22ns1*Nd(ksi(I),2);
    g22ns1(3)=g22ns1(3)+ig22ns1*Nd(ksi(I),3);
        % ns2 com qsi1
    hns2(1)=hns2(1)+ihns2*Nd(qsi1,1);
    hns2(2)=hns2(2)+ihns2*Nd(qsi1,2);
    hns2(3)=hns2(3)+ihns2*Nd(qsi1,3);

    g11ns2(1)=g11ns2(1)+ig11ns2*Nd(qsi1,1);
    g11ns2(2)=g11ns2(2)+ig11ns2*Nd(qsi1,2);
    g11ns2(3)=g11ns2(3)+ig11ns2*Nd(qsi1,3);
    g22ns2(1)=g22ns2(1)+ig22ns2*Nd(qsi1,1);
    g22ns2(2)=g22ns2(2)+ig22ns2*Nd(qsi1,2);
    g22ns2(3)=g22ns2(3)+ig22ns2*Nd(qsi1,3);
        % ns3 com qsi2
    hns3(1)=hns3(1)+ihns3*Nd(qsi2,1);
    hns3(2)=hns3(2)+ihns3*Nd(qsi2,2);
    hns3(3)=hns3(3)+ihns3*Nd(qsi2,3);

    g11ns3(1)=g11ns3(1)+ig11ns3*Nd(qsi2,1);
    g11ns3(2)=g11ns3(2)+ig11ns3*Nd(qsi2,2);
    g11ns3(3)=g11ns3(3)+ig11ns3*Nd(qsi2,3);
    g22ns3(1)=g22ns3(1)+ig22ns3*Nd(qsi2,1);
    g22ns3(2)=g22ns3(2)+ig22ns3*Nd(qsi2,2);
    g22ns3(3)=g22ns3(3)+ig22ns3*Nd(qsi2,3);

    % Integracao singular-singular (logaritma): Mn, w e d2wdndni
        % s1 com qsi1l
    hs1(1)=hs1(1)+ihs1*Nd(qsi1l,1);
    hs1(2)=hs1(2)+ihs1*Nd(qsi1l,2);
    hs1(3)=hs1(3)+ihs1*Nd(qsi1l,3);        

    g11s1(1)=g11s1(1)+ig11s1*Nd(qsi1l,1);
    g11s1(2)=g11s1(2)+ig11s1*Nd(qsi1l,2);
    g11s1(3)=g11s1(3)+ig11s1*Nd(qsi1l,3);
    g22s1(1)=g22s1(1)+ig22s1*Nd(qsi1l,1);
    g22s1(2)=g22s1(2)+ig22s1*Nd(qsi1l,2);
    g22s1(3)=g22s1(3)+ig22s1*Nd(qsi1l,3);
        % s2 com qsi2l
    hs2(1)=hs2(1)+ihs2*Nd(qsi2l,1);
    hs2(2)=hs2(2)+ihs2*Nd(qsi2l,2);
    hs2(3)=hs2(3)+ihs2*Nd(qsi2l,3);

    g11s2(1)=g11s2(1)+ig11s2*Nd(qsi2l,1);
    g11s2(2)=g11s2(2)+ig11s2*Nd(qsi2l,2);
    g11s2(3)=g11s2(3)+ig11s2*Nd(qsi2l,3);
    g22s2(1)=g22s2(1)+ig22s2*Nd(qsi2l,1);
    g22s2(2)=g22s2(2)+ig22s2*Nd(qsi2l,2);
    g22s2(3)=g22s2(3)+ig22s2*Nd(qsi2l,3);

    % Parte regular: Mn e d2wdndni (w nao tem parte reg)      
    hreg(1)=hreg(1)+ihreg*Nd(ksi(I),1);
    hreg(2)=hreg(2)+ihreg*Nd(ksi(I),2);
    hreg(3)=hreg(3)+ihreg*Nd(ksi(I),3);

    g22reg(1)=g22reg(1)+ig22reg*Nd(ksi(I),1);
    g22reg(2)=g22reg(2)+ig22reg*Nd(ksi(I),2);
    g22reg(3)=g22reg(3)+ig22reg*Nd(ksi(I),3);
end

% Sao zeros de fatos % Atribuicao na rotina externa n2
% h(1,1)=0; % Vn
% h(2,2)=0; % dMndni
% g(1,2)=0; % dwdn
% g(2,1)=0; % dwdni

% Mn Singularidade Logaritmica h(1,2)
hMn(1)=hreg(1)+hns1(1)+hns2(1)+hns3(1)+hs1(1)+hs2(1);
hMn(2)=hreg(2)+hns1(2)+hns2(2)+hns3(2)+hs1(2)+hs2(2);
hMn(3)=hreg(3)+hns1(3)+hns2(3)+hns3(3)+hs1(3)+hs2(3);

% w Singularidade Logaritmica g(1,1)
g11(1)=g11reg(1)+g11ns1(1)+g11ns2(1)+g11ns3(1)+g11s1(1)+g11s2(1);
g11(2)=g11reg(2)+g11ns1(2)+g11ns2(2)+g11ns3(2)+g11s1(2)+g11s2(2);
g11(3)=g11reg(3)+g11ns1(3)+g11ns2(3)+g11ns3(3)+g11s1(3)+g11s2(3);

% d2wdndni Singularidade Logaritmica g(2,2)
g22(1)=g22reg(1)+g22ns1(1)+g22ns2(1)+g22ns3(1)+g22s1(1)+g22s2(1);
g22(2)=g22reg(2)+g22ns1(2)+g22ns2(2)+g22ns3(2)+g22s1(2)+g22s2(2);
g22(3)=g22reg(3)+g22ns1(3)+g22ns2(3)+g22ns3(3)+g22s1(3)+g22s2(3);

% dVndni HipersingularHadamard dVndni=(1+ni)/(4*pi*R^2);
ep=2e-16; % tolerancia de comparacao
if x2>(ep+(x1+x3)/2)||x2<(-ep+(x1+x3)/2)||y2>(ep+(y1+y3)/2)||y2<(-ep+(y1+y3)/2)
    disp('ERROn3_h1g1Sing12:Foi encontrado um elemento curvo')
    x2-(x1+x3)/2
    y2-(y1-y3)/2
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
