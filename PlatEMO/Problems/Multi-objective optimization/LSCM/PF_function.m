function [Obj] = PF_function(PopDec,G,PF_form,Q_form,N,M)

if PF_form == 1
    if Q_form == 1
        Obj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
    elseif Q_form == 2
        Obj = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
    elseif Q_form == 3
        G = sum(G,2);
        Obj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
    end
    
elseif PF_form == 2
    if Q_form == 1
        Obj =  (1+G).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
    elseif Q_form == 2
        Obj =  (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
    elseif Q_form == 3
        G = sum(G,2);
        Obj =  (1+G).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
    end
elseif PF_form == 3
    if Q_form == 1
        Obj(:,1:M-1) = PopDec(:,1:M-1);
        G = 1 + G(:,end);
        Obj(:,M)     = (1+G).*(M-sum(Obj(:,1:M-1)./(1+repmat(G,1,M-1)).*(1+sin(3*pi.*Obj(:,1:M-1))),2));
    end
end
end

