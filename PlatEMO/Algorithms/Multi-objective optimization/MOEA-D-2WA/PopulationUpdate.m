function Population = PopulationUpdate(Population, N, initialE, epsn, Zmin)
% Update the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    nCon   = size(Population.cons, 2);
    ConVio = max(0,Population.cons);
    if sum(sum(ConVio<=epsn, 2)==nCon) >= N
        %% Selection among epsilon-feasible solutions
        tmp        = sum(ConVio<=epsn, 2)==nCon;
        Population = Population(1:end, tmp);
        PopCv      = sum(max(0, Population.cons)./initialE, 2)./nCon;
        
        %% calculate angle between each two solutions objective values
        Next(1:size(Population, 2)) = true;
        PopObj                      = Population(Next).objs;
        Delete                      = LastSelection(PopObj, sum(Next) - N, Zmin, PopCv);
        Temp                        = find(Next);
        Next(Temp(Delete))          = false;
        Population                  = Population(Next);
    else
        %% Selection including epsilon-infeasible solutions
        CV         = sum(max(0, Population.cons)./initialE, 2)./nCon;
        [~, rank]  = sort(CV);
        % Population for next generation
        Population = Population(rank(1:N));
    end
end

function Delete = LastSelection(PopObj, K, Zmin, PopCv)
    %% calculate angle between each two solutions based on objective values
    N      = size(PopObj, 1);
    PopObj = (PopObj - repmat(Zmin, N, 1))./(repmat(max(PopObj), N, 1) - repmat(Zmin, N, 1));
    Cosine = 1 - pdist2(PopObj, PopObj, 'cosine');
    Cosine = Cosine.*(1 - eye(size(PopObj, 1)));

    % Environmental selection
    Delete  = false(1, N);
    % Delete K solutions one by one with smallest angle
    while sum(Delete) < K
        [Jmin_row, Jmin_column] = find(Cosine==max(max(Cosine)));
        j      = randi(length(Jmin_row));
        Temp_1 = Jmin_row(j);
        Temp_2 = Jmin_column(j);
        if PopCv(Temp_1) > PopCv(Temp_2)
            Delete(Temp_1)    = true;
            Cosine(:, Temp_1) = 0;
            Cosine(Temp_1, :) = 0;
        else
            Delete(Temp_2)    = true;
            Cosine(:, Temp_2) = 0;
            Cosine(Temp_2, :) = 0;
        end 
    end
end