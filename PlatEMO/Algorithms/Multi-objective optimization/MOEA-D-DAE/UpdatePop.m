function [Population,pi] = UpdatePop(Population,P,Offspring,epsilon,z,W,sigma,pi,avg_fit)
% Update population 

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
        g_pi  = max(repmat(abs(Population(i).objs-z),1,1).*W(i,:),[],2);
        if (sum(max(Population(i).cons,0),2)<=epsilon) && (sum(max(Offspring.cons,0),2)<epsilon)
            fitnesso  = g_o;
            fitnesspi = g_pi;
        else
            fitnesso  = sigma*g_o+(1-sigma)*avg_fit*sum(max(Offspring.cons,0),2);
            fitnesspi = sigma*g_pi+(1-sigma)*avg_fit*sum(max(Population(i).cons,0),2);
        end
        if fitnesso < fitnesspi
            Population(i) = Offspring;
            delta = (fitnesspi-fitnesso)/fitnesspi;
            if delta > 0.001
                pi(i) = 1;
            else
                pi(i) = 0.95+(0.05*delta/0.001)*pi(i); 
            end
            c = c+1;
        end
        P(index(1)) = [];
    end
end