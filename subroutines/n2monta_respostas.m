function [w,Vn,dwdn,Mn,wc,Rc]=n2monta_respostas(x,CDC,CANTOS)

for I = 1:size(CDC,1)  % Percorre os elementos
    
    tipoCDC=CDC(I,2); % Tipo da condicao de contorno no elemento I
    
    if tipoCDC==1 % elemento engastado
        w([3*I-2,3*I-1,3*I])=CDC(I,3:2:7);
        dwdn([3*I-2,3*I-1,3*I])=CDC(I,4:2:8);
        Vn([3*I-2,3*I-1,3*I])=x([6*I-5,6*I-3,6*I-1]);
        Mn([3*I-2,3*I-1,3*I])=x([6*I-4,6*I-2,6*I-0]);
        
    elseif tipoCDC==2   % elemento apoiado
        w([3*I-2,3*I-1,3*I])=CDC(I,3:2:7);
        Mn([3*I-2,3*I-1,3*I])=CDC(I,4:2:8);
        Vn([3*I-2,3*I-1,3*I])=x([6*I-5,6*I-3,6*I-1]);
        dwdn([3*I-2,3*I-1,3*I])=x([6*I-4,6*I-2,6*I-0]);
        
    elseif tipoCDC==3 % elemento guia
        Vn([3*I-2,3*I-1,3*I])=CDC(I,3:2:7);
        dwdn([3*I-2,3*I-1,3*I])=CDC(I,4:2:8);
        w([3*I-2,3*I-1,3*I])=x([6*I-5,6*I-3,6*I-1]);
        Mn([3*I-2,3*I-1,3*I])=x([6*I-4,6*I-2,6*I-0]);
        
    elseif tipoCDC==4 % elemento livre
        Vn([3*I-2,3*I-1,3*I])=CDC(I,3:2:7);
        Mn([3*I-2,3*I-1,3*I])=CDC(I,4:2:8);
        w([3*I-2,3*I-1,3*I])=x([6*I-5,6*I-3,6*I-1]);
        dwdn([3*I-2,3*I-1,3*I])=x([6*I-4,6*I-2,6*I-0]);
        
    else disp('ERRO:tipoCDC deve estar entre 1 e 4')
    end
end

for I = 1:size(CANTOS,1)
    J=I+6*size(CDC,1); % indica a parte final do vetor x
    if CANTOS(I,5)==1 %preso
        wc(I)=CANTOS(I,6);
        Rc(I)=x(J);
    elseif CANTOS(I,5)==2 %solto
        Rc(I)=CANTOS(I,6);
        wc(I)=x(J);
    else disp('ERRO:tipoCDCcanto deve ser 1 ou 2')
    end
end