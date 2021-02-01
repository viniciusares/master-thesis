function [w]=timoEALA(x,y,a,b,final)

% valor de "final" ate 11 vai bem

q=1;
D=1;
%a=1; b=1;
v=0.3;

soma1=0;
w2=0;
for m=1:2:final
soma1=soma1+1/m^5*sin(m*pi*x/a);


betam=m*pi*b/a;
Am=-4/(pi^5*m^5);
Bm=4/(pi^5*m^5)*((3+v)*(1-v)*(cosh(betam))^2+2*v*cosh(betam)-...
    v*(1-v)*betam*sinh(betam)-(1-v^2))/...
    ((3+v)*(1-v)*(cosh(betam))^2+(1-v)^2*betam^2+(1+v)^2);
Cm=4/(pi^5*m^5)*((3+v)*(1-v)*sinh(betam)*cosh(betam)+v*(1+v)*sinh(betam)-...
    v*(1-v)*betam*cosh(betam)-(1-v)^2*betam)/...
    ((3+v)*(1-v)*(cosh(betam))^2+(1-v)^2*betam^2+(1+v)^2);
Dm=-Cm;
Ym=q*a^4/D*(Am*cosh(m*pi*y/a)+Bm*m*pi*y/a*sinh(m*pi*y/a)+...
    Cm*sinh(m*pi*y/a)+Dm*m*pi*y/a*cosh(m*pi*y/a));
w2=w2+Ym*sin(m*pi*x/a);
end


w1=4*q*a^4/(pi^5*D)*soma1;

w=w1+w2;
