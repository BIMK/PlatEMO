function tA = UpdatetA(tA,P,Offspring,z,W,avg_fit)
% Update temparory archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    c  = 0;
    nr = 2; 
    while c~=nr && size(P,2)~=0
        index = randperm(size(P,2));
        i     = P(index(1));
        g_o   = max(repmat(abs(Offspring.objs-z),1,1).*W(i,:),[],2);
        g_pi  = max(repmat(abs(tA(i).objs-z),1,1).*W(i,:),[],2);
        e     = 1/size(tA,2);
        if (e*g_o+(1-e)*avg_fit*sum(max(Offspring.cons,0),2))<(e*g_pi+(1-e)*avg_fit*sum(max(tA(i).cons,0),2))
            tA(i) = Offspring;
        end
        c = c+1;
        P(index(1)) = [];
    end
end