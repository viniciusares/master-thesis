function [H1int,G1int,R1int,C1int,u,q]=...
    wIntern(P_int,ELGEO,GEO,CANTOS,w,dwdn,Vn,Mn,wc,Rc,D,v)

% Calcula a deflexao nos pontos internos

n_pint=size(P_int,1); % Numero de pontos internos
n_elem=size(ELGEO,1); % Numero de elementos
n_cantos=size(CANTOS,1); % Numero de cantos
H1int=zeros(n_pint,6*n_elem); % Inicializa matriz H1_int
G1int=zeros(n_pint,6*n_elem); % Inicializa matriz G1_int
R1int=zeros(n_pint,n_cantos); % Inicializa matriz R1_int
C1int=zeros(n_pint,n_cantos); % Inicializa matriz C1_int
u=zeros(1,2*size(w,2));
q=zeros(1,2*size(Vn,2));

for I=1:n_pint % Loop sobre os pontos internos
    
    X_fonte=P_int(I,2:3); % Coordenadas do ponto fonte
    
    for J=1:n_elem
        % Coordenadas do elemento
        X1=GEO(ELGEO(J,2),2:3);
        X2=GEO(ELGEO(J,3),2:3); 
        X3=GEO(ELGEO(J,4),2:3);
        
        [h,g]=n3_h2g2Regular(X_fonte,X1,X2,X3,D,v);
        
        H1int(I,6*J-5:6*J)=h(:); % sem eq hipersing
        G1int(I,6*J-5:6*J)=g(:); % sem eq hipersing
    end
    
    for J = 1 : n_cantos % Percorre os cantos, R1C1
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
        
        [r1,c1] = n3_r2c2Regular(X_fonte,Xcj,Nj,NJ,Sj,SJ,D,v);
        
        % Organiza componentes vetor por vetor (2x1)
        R1int(I,J)=r1; %contem Rc*
        C1int(I,J)=c1; %contem wc*
    end
end

u(1:2:end)=w(:);    % Colunas impares
u(2:2:end)=dwdn(:); % Colubas pares
for I=1:n_cantos
    J=I+6*size(ELGEO,1);
    u(J)=wc(I);
end

q(1:2:end)=Vn(:);   % Colunas impares
q(2:2:end)=Mn(:);   % Colunas pares
for I=1:n_cantos
    J=I+6*size(ELGEO,1);
    q(J)=Rc(I);
end
