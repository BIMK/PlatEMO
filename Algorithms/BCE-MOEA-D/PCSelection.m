function [PC,nND] = PCSelection(PC,N)
% PC selection and population maintenance in BCE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% PC selection
    PC    = PC(NDSort(PC.objs,1)==1);
    PC    = PC(randperm(length(PC)));
    PCObj = PC.objs;
    nND   = length(PC);
    
    %% Population maintenance
    if length(PC) > N
        % Normalization
        fmax  = max(PCObj,[],1);
        fmin  = min(PCObj,[],1);
        PCObj = (PCObj-repmat(fmin,nND,1))./repmat(fmax-fmin,nND,1);
        % Determine the radius of the niche
        d  = pdist2(PCObj,PCObj);
        d(logical(eye(length(d)))) = inf;
        sd = sort(d,2);
        r  = mean(sd(:,min(3,size(sd,2))));
        R  = min(d./r,1);
        % Delete solution one by one
        while length(PC) > N
            [~,worst]  = max(1-prod(R,2));
            PC(worst)  = [];
            R(worst,:) = [];
            R(:,worst) = [];
        end
    end
end