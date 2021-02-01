switch cruzougrade
    case 1 % So serve para EALA
        % 3*3*11-(3*((11-1)/2)+1) = 99-16=83
        nohcentral=3*3*n-(3*((n-1)/2)+1); % Para descobrir qual noh central
        (w_dwdn(nohcentral)*D-timoEALA(0.5,1,11))/timoEALA(0.5,1,11)
        (w_int(5)*D-timoEALA(0.5,0.5,11))/timoEALA(0.5,0.5,11) % sempre 5
    case 2
        figure (1)
        xpasso=(xmax-xmin)/(xnpts-1); ypasso=(ymax-ymin)/(ynpts-1);
        xgrid=xmin:xpasso:xmax; ygrid=ymin:ypasso:ymax;
        [xx,yy]=meshgrid(xgrid,ygrid);
        Wgrid=FgeraWgrid(w_int,xnpts,ynpts);
        subplot(1,2,1)
        surf(xx,yy,Wgrid*D/(cargcel*LadoX^2*LadoY^2))
        title('ProgKirch w_{adi}=w*D/(cargcel*LadoX^2*LadoY^2)');
        xlabel('eixoX'); ylabel('eixoY')
        foo = get(gca,'dataaspectratio');
        set(gca,'dataaspectratio',[foo(1) foo(1) 0.2*foo(1)]);
        
        for J=1:size(xx,2)
            for I=1:size(yy,1)
                WgridTimo(I,J)=timoEALA(xx(I,J),yy(I,J),LadoX,LadoY,7);
            end
        end
        subplot(1,2,2)
        surf(xx,yy,WgridTimo)
        title('WtimoEALA'); xlabel('eixoX'); ylabel('eixoY')
        foo = get(gca,'dataaspectratio');
        set(gca,'dataaspectratio',[foo(1) foo(1) 1*foo(1)]);
%================================================================
%         figure (2) % Perfil de deflexoes EALA Assimetrico
%         Jxcte=xnpts/2+0.5; % x=0.5 Sobre o Eixo de simetria
%         deflsFy=Wgrid(:,Jxcte);
%         adiFy=D*deflsFy/(cargcel*LadoX^2*LadoY^2);
%         subplot(2,1,1)
%         plot(ygrid,adiFy,'r')
%         title('w_{adi}(y): do E ate L. x=0.5(cte)')
%         hold on
%         deflsFyTIMO=zeros(size(ygrid));
%         for I=1:length(ygrid)
%             deflsFyTIMO(I)=...
%                 FaaaTimoInv(0.5*LadoX,ygrid(I),LadoX,LadoY,h,E,cargcel,11);
%         end
%         adiFyT=D*deflsFyTIMO/(cargcel*LadoX^2*LadoY^2);
%         plot(ygrid,adiFyT,'b');
%         legend('matlab','timoshenko')
% 
%         % Perfil de deflexoes EALA Simetrico
%         Iycte=ynpts/2+0.5; % y=0.5 Meio-caminho entre E e L
%         % Cruza o eixo de simetria, de A ate A
%         deflsFx=Wgrid(Iycte,:);
%         adiFx=D*deflsFx/(cargcel*LadoX^2*LadoY^2);
%         subplot(2,1,2)
%         plot(xgrid,adiFx','r')
%         title('w_{adi}(x): de E ate E, simetrico. y=0.5(cte)')
%         hold on
%         deflsFxTIMO=zeros(size(xgrid));
%         for I=1:length(xgrid)
%             deflsFxTIMO(I)=...
%                 FaaaTimoInv(xgrid(I),0.5*LadoY,LadoX,LadoY,h,E,cargcel,11);
%         end
%         adiFxT=D*deflsFxTIMO/(cargcel*LadoX^2*LadoY^2);
%         plot(xgrid,adiFxT,'b');
%         legend('matlab','timoshenko')
%         %==============================================
%         figure (3)
%         subplot(1,2,1)
%         plot(ygrid,adiFy'./adiFyT);
%         title(['Erro de w(y). w(y)/w_{timo}(y)'])
%         
%         subplot(1,2,2)
%         plot(xgrid,adiFx./adiFxT);
%         title(['Erro de w(x). w(x)/w_{timo}(x)'])
    case 3 % Cruz AAAA
        w_vini=w_int(1:5).*D;
        w_timo_w_vini=[w_timo,w_vini]
        erroAAAA=(w_vini-w_timo)./w_timo
end