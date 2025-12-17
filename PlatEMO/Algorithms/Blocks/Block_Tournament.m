classdef Block_Tournament < BLOCK
% Tournament selection
% nParents --- 100 --- Number of parents generated
% upper    ---   2 --- Max number of k for k-tournament

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nParents;       % <hyperparameter> Number of generated parents 
        nTournament;    % <parameter> Number of candidate solutions for selecting a parent
    end
    methods
        %% Default settings of the block
        function obj = Block_Tournament(nParents,upper)
            obj.nParents = nParents;    % Hyperparameter
            obj.lower    = 1;           % Lower bounds of parameters
            obj.upper    = upper;      	% Upper bounds of parameters
            % Randomly set the parameters
            obj.parameter = unifrnd(obj.lower,obj.upper);
            obj.ParameterAssign();
        end
        %% Assign parameters to variables
        function ParameterAssign(obj)
            obj.nTournament = round(obj.parameter(1));
        end
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            Population = obj.Gather(Problem,Precursors,Ratio,1,1);
            if obj.nTournament == 1
                obj.output = Population(randi(end,1,obj.nParents));
            elseif Problem.M == 1   % For single-objective optimization
                obj.output = Population(TournamentSelection(obj.nTournament,obj.nParents,FitnessSingle(Population)));
            else                    % For multi- and many-objective optimization
                [FrontNo,Dis] = CalFitness(Population);
                obj.output    = Population(TournamentSelection(obj.nTournament,obj.nParents,FrontNo,Dis));
            end
        end
    end
end

function [FrontNo,Dis] = CalFitness(Population)
    Dis = zeros(length(Population),1);
    [FrontNo,maxF] = NDSort(Population.objs,Population.cons,inf);
    for front = 1 : maxF
        PopObj = Population(FrontNo==front).objs;
        if size(PopObj,2) < 4
            Distance = pdist2(PopObj,PopObj);
            Distance(logical(eye(size(PopObj,1)))) = inf;
        else
            Distance = inf(size(PopObj,1));
            for i = 1 : size(PopObj,1)
                SPopObj = max(PopObj,repmat(PopObj(i,:),size(PopObj,1),1));
                for j = [1:i-1,i+1:size(PopObjAll,1)]
                    Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
                end
            end
        end
        Distance = sort(Distance,2);
        Dis(FrontNo==front) = sum(1./Distance(:,1:min(end,2)),2);
    end
end