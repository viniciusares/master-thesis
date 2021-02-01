clear all
close all
clc

maxtimo=2001;
mmax=51; %51 da 6 algarismos significativos de precisao
nmax=51; %acima de 1001 fica lento
a=1;
b=1;
q0=1;
D=1;
x=1;
y=1;

w=0;
if mmax<=maxtimo && nmax<=maxtimo
    for m=1:2:mmax
        for n=1:2:nmax

            w=w+sin(m*pi*x/a)*sin(n*pi*y/b)/(m*n*(m^2/a^2+n^2/b^2)^2);

        end
    end
else
    disp('Error:maxtimo exceeded')
end

format long
w=w*16*q0/(pi^6*D)