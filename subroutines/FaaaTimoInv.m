function w=FaaaTimoInv(x,y,a,b,h,E,q0,termos)

maxtimo=2001;
mmax=termos; nmax=termos; %51 da 6 algarismos de precisao, >1001 fica lento
% a=1; b=1;
v=0.3;
D=E*h^3/(12*(1-v^2));
%x=0.1; y=0.5;

w=0;
if mmax<=maxtimo && nmax<=maxtimo
    
    for I=1:2:mmax
        m=mmax+1-I; % soma primeiro os menores termos
        
        for J=1:2:nmax
            n=nmax+1-J; % soma primeiro os menores termos

            w=w+sin(m*pi*x/a)*sin(n*pi*y/b)/(m*n*(m^2/a^2+n^2/b^2)^2);

        end
    end
    
else
    disp('Error:maxtimo exceeded')
end

format long
w=w*16*q0/(pi^6*D);