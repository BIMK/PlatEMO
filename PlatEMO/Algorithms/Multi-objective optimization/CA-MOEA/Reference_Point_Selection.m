function  [reference_population]=Reference_Point_Selection(MaxFnorm,last_population,Ref,K,M)
% Clustering-based environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yicun Hua

    last_Num = size(MaxFnorm,1);
    Ref_Num = size(Ref,1);
    dis = pdist2(MaxFnorm,Ref);
    
    for i = 1:last_Num
        [~,dmin_p] = min(dis(i,:),[],2);
        
        % Individual i assigned to Ref(dmin_p,:)
        MaxFnorm(i,(M+1)) = dmin_p; 
        
        % Distance between individual i and Ref(dmin_p,:)
        MaxFnorm(i,(M+2)) = dis(i,dmin_p);
    end
    
    % Individuals are labeled 3(to distinguish from 1 and 2)
    MaxFnorm(:,(M+3)) = 3;
    for i = 1:Ref_Num
        if sum(MaxFnorm(:,(M+1)) == i) ~= 0
            a = find(MaxFnorm(:,(M+1)) == i);
            a = a';
            [~,b] = sort(MaxFnorm(a,(M+2)));
            
            % Individuals closest to reference points are labeled 1
            MaxFnorm(a(b(1)),(M+3)) = 1;   
            
            % Individuals in crowd area are labeled 2
            if sum(MaxFnorm(:,(M+1)) == i) > 3
                a(b(1)) = [];
                MaxFnorm(a,(M+3)) = 2; 
            end
        end
    end
    
    t1 = sum(MaxFnorm(:,(M+3)) == 1);
    t2 = sum(MaxFnorm(:,(M+3)) == 2); 
    
    if t1 <= K
        % If individuals labeled 1 are not enough, we choose some 2, then 3
        d = find(MaxFnorm(:,(M+3))==1);
        d = d';
        reference_population(1:t1) = d;
        if t2 >= (K-t1)
            % If individuals labeled 2 are enough
            d = find(MaxFnorm(:,(M+3))==2);
            d = d';
            y = randperm(length(d));
            reference_population((t1+1):K) = d(y(1:(K-t1)));
        else
            % If individuals labeled 2 are not enough, we choose some 3
            d = find(MaxFnorm(:,(M+3)) == 2);
            d = d';
            reference_population((t1+1):(t1+t2)) = d;
            d = find(MaxFnorm(:,(M+3)) == 3);
            d = d';
            y = randperm(length(d));
            reference_population((t1+t2+1):K) = d(y(1:(K-t1-t2)));
        end
    else
        % If individuals labeled 1 are enough
        d = find(MaxFnorm(:,(M+3)) == 1);
        d = d';
        y = randperm(length(d));
        reference_population(1:K) = d(y(1:K));
    end
    reference_population = last_population(reference_population);
end