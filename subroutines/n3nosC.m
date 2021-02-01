function [nosC]=n3nosC(x1,y1,x2,y2,x3,y3)

nosC=zeros(3,2);
for I=1:3
    ksi=-4/3+I*2/3;
    phi1c=1/2*ksi*(ksi-1);  %funcao de forma 1
    phi2c=1-ksi^2;          %funcao de forma 2
    phi3c=1/2*ksi*(ksi+1);  %funcao de forma 3
    nosC(I,1)=phi1c*x1+phi2c*x2+phi3c*x3;
    nosC(I,2)=phi1c*y1+phi2c*y2+phi3c*y3;
end