function [A,B,C,D,AmC,BmD]=testetimoEALA(x,y,m)
% teste timoEALA

q=1;
D=1;
a=1; b=1;
v=0.3;

% Calculo de w1
soma1=1/m^5*sin(m*pi*x/a);
w1=4*q*a^4/(pi^5*D)*soma1;

% Calculo de w2
betam=m*pi*b/a;
Am=-4/(pi^5*m^5);
Bm=4/(pi^5*m^5)*((3+v)*(1-v)*(cosh(betam))^2+2*v*cosh(betam)-...
    v*(1-v)*betam*sinh(betam)-(1-v^2))/...
    ((3+v)*(1-v)*(cosh(betam))^2+(1-v)^2*betam^2+(1+v)^2);
Cm=4/(pi^5*m^5)*((3+v)*(1-v)*sinh(betam)*cosh(betam)+v*(1+v)*sinh(betam)-...
    v*(1-v)*betam*cosh(betam)-(1-v)^2*betam)/...
    ((3+v)*(1-v)*(cosh(betam))^2+(1-v)^2*betam^2+(1+v)^2);
Dm=-Cm;
Ymch=cosh(m*pi*y/a);
Ymsh=sinh(m*pi*y/a);
A=Am*Ymch;
B=Bm*m*pi*y/a*Ymsh;
C=Cm*Ymsh;
D=Dm*m*pi*y/a*Ymch;
Ym=q*a^4/D*(A+B+C+D);

w2seno=sin(m*pi*x/a);
w2=Ym*w2seno;

w=w1+w2;

AmC=A+C;
BmD=B+D;
