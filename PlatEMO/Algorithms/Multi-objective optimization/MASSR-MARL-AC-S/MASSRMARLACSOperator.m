function [OffDec,OffMask] = MASSRMARLACSOperator(Problem,ParentDec,ParentMask,Fitness,Mask,operator,Memory,pv,fv,progress)
%MASSRMARLACSOperator - Generate latent decisions and sparse masks.
%
%   operator 1 uses GLP-GA and operator 2 uses MSKEA. Parent decisions
%   remain separate from their binary masks; the caller evaluates
%   OffDec.*OffMask.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N                              = size(ParentDec,1)/2;
    D                              = Problem.D;
    first                          = ParentMask(1:N,:);
    second                         = ParentMask(N+1:end,:);
    [~,order]                      = sort(Fitness);
    group                          = ones(1,D);
    group(order(floor(D/2)+1:end)) = 2;

    % GLP-GA: grouped decision variation and mask updates
    if operator == 1
        OffDec      = MASSRMARLACSGLP(Problem,ParentDec(1:N,:),ParentDec(N+1:end,:),4);
        OffMask     = first;
        diversity   = max(1,numel(unique(Memory(max(1,end-9):end))));
        probability = min(1,100*mean(Mask(:))/(D*diversity^2));
        for i = 1 : N
            selected           = randi(2);
            changed            = xor(first(i,:),second(i,:)) & group == selected;
            OffMask(i,changed) = rand < 0.5;
            changed            = rand(1,D) < probability;
            OffMask(i,changed) = rand < 0.5;
        end
    else
        % MSKEA: decision crossover and score-guided mask changes
        if any(Problem.encoding~=4)
            OffDec = OperatorGAhalf(Problem,ParentDec);
        else
            OffDec = ones(N,D);
        end
        OffMask = first;
        for i = 1 : N
            scale  = 0.8 + 0.2*(group==randi(2));
            score  = pv.*scale;
            spread = fv.*scale;
            if rand < 0.5
                indices = find(first(i,:)&~second(i,:));
                if ~isempty(indices)
                    OffMask(i,indices(TournamentSelection(2,1,-score(indices)))) = false;
                end
            else
                indices = find(~first(i,:)&second(i,:));
                if ~isempty(indices)
                    OffMask(i,indices(TournamentSelection(2,1,score(indices)))) = true;
                end
            end
            if rand < 1 - progress
                indices = find(OffMask(i,:)~=(spread>0));
                sign    = rand < 0.5;
                if sign
                    rank = -spread;
                else
                    rank = spread;
                end
            else
                sign = rand >= 0.5;
                if sign
                    indices = find(~OffMask(i,:));
                    rank    = score;
                else
                    indices = find(OffMask(i,:));
                    rank    = -score;
                end
            end
            if ~isempty(indices)
                OffMask(i,indices(TournamentSelection(2,1,rank(indices)))) = sign;
            end
        end
    end
    OffDec(:,Problem.encoding==4) = 1;
    OffMask                       = logical(OffMask);
end