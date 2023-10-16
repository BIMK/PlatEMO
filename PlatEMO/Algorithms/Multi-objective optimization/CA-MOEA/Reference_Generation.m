function [Ref] = Reference_Generation(MaxFnorm, M, K)
% Clustering-based reference points generation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yicun Hua

    % Calculate boundary individuals: extreme
    [~,maxobj] = max(MaxFnorm,[],1);
    [~,minobj] = min(MaxFnorm,[],1);
    ex = unique([maxobj,minobj]);
    extreme = unique(MaxFnorm(ex,:),'rows');
    
    E = size(extreme,1);    
    Nc = K-E;
    
    if Nc>0
        % Clustering
        T = clusterdata(MaxFnorm,'maxclust',Nc,'distance','euclidean','linkage','ward');
        
        % Calculate cluster centers
        for i = 1:Nc

            p = find(T == i);
            pn = length(p);
 
            Ref(i,1:M) = sum(MaxFnorm(p,:),1)/pn;

            Ref(i,M+1) = i;   
            
        end
        
        Ref = unique(Ref,'rows');
        rpr = size(Ref,1);
        
        % If cluster centers are not enough, 
        % randomly select some individuals as reference points
        if rpr < Nc
            Ncrp = Nc-rpr;
            Ref(rpr+1:Nc,:) = MaxFnorm(randperm(size(MaxFnorm,1),Ncrp),:);
        end
       
        Ref(:,M+1) = [];
        
        % Use boundary individuals as reference points
        for er = 1:E
            Ref(Nc+er,1:M) = extreme(er,1:M);
        end
    else
        % Use boundary individuals as reference points
        for er = 1:K
            Ref(er,:) = extreme(er,:);
        end
    end
end