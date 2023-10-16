function Archive = ArchiveUpdate(Archive,N)
% Archive Update

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if isempty(Archive)
        return;
    else       
        nND = length(Archive);
        while length(Archive) > N
            % Normalization
            PCObj = Archive.objs;
            fmax  = max(PCObj,[],1);
            fmin  = min(PCObj,[],1);
            PCObj = (PCObj-repmat(fmin,nND,1))./repmat(fmax-fmin,nND,1);
            % Determine the radius of the niche
            d  = pdist2(PCObj,PCObj);
            d(logical(eye(length(d)))) = inf;
            sd = sort(d,2);
            r  = median(sd(:,min(size(PCObj,2),size(sd,2))));
            R  = min(d./r,1);
            % Delete solution one by one
            while length(Archive) > N
                [~,worst]  = max(1-prod(R,2));
                Archive(worst)  = [];
                R(worst,:) = [];
                R(:,worst) = [];
            end         
        end
    end
end