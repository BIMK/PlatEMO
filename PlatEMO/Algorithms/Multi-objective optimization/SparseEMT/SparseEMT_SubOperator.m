function [OffDec,OffMask]  = SparseEMT_SubOperator(Problem,ParentDec,ParentMask,form,Dim)
%% Crossover and mutation operators of SparseEA

    % Parameter setting
    [N,D]       = size(ParentDec);
    
    % Crossover and mutation
    if N > 1
        if form == 3
            Parent1Mask = ParentMask(1:floor(N/2),:);
            Parent2Mask = ParentMask(floor(N/2)+1:end,:);

            % Crossover for mask
            k = rand(N/2,D) < 0.5;
            k(repmat(rand(N/2,1)>1,1,D)) = false;
            OffMask    = Parent1Mask;
            OffMask(k) = Parent2Mask(k);

            % Mutation for mask
            for i = 1 : floor(N/2)
                Site = rand(1,D) < 1/D;
                OffMask(i,Site) = ~OffMask(i,Site);
            end

            % The same position of the parent mask is inherited
            for i = 1 : D
                if Parent1Mask(:,i) == Parent2Mask(:,i)
                    OffMask(:,i) = Parent1Mask(:,i);
                end
            end

        else
            OffMask = ones(floor(N/2),D);
        end

    else
        % Mutation for mask
        OffMask = ParentMask;
        if form == 3
            Site = rand(1,D) < 1/D;
            OffMask(:,Site) = ~OffMask(:,Site);
        end
    end



    
    % Crossover and mutation for dec
    if N > 1
        if any(Problem.encoding~=4)
            OffDec = SparseEMT_OperatorGAhalf(Problem,ParentDec,Dim);
            OffDec(:,Problem.encoding==4) = 1;
        else
            OffDec = ones(floor(N/2),D);
        end
    else
        if any(Problem.encoding~=4)
            OffDec = OperatorMutate(Problem,ParentDec,Dim);
            OffDec(:,Problem.encoding==4) = 1;
        else
            OffDec = ParentDec;
        end
    end
end


function Offspring = OperatorMutate(Problem,Parent,Dim,Parameter)
%% Mutation for real and integer variables

    if nargin > 3
        [proM,disM] = deal(Parameter{:});
    else
        [proM,disM] = deal(1,20);
    end
    D     = length(Dim);
    Lower = repmat(Problem.lower(Dim),1,1);
    Upper = repmat(Problem.upper(Dim),1,1);
    Site  = rand(1,D) < proM/D;
    mu    = rand(1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Parent,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end