classdef Block_Crossover < BLOCK
% Unified crossover for real variables
% nParents --- 2 --- Number of parents generating one offspring
% nSets    --- 5 --- Number of parameter sets

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nParents;	% <hyperparameter> Number of parents for generating an offspring
        nSets;      % <hyperparameter> Number of weight sets
        Weight;     % <parameter> Weight sets
        Fit;        % <parameter> Probability of using each weight set
    end
    methods
        %% Default settings of the block
        function obj = Block_Crossover(nParents,nSets)
            obj.nParents = nParents;	% Hyperparameter
            obj.nSets    = nSets;     	% Hyperparameter
            obj.lower    = repmat([0 -1 1e-20],1,(nParents-1)*nSets);	% Lower bounds of parameters
            obj.upper    = ones(1,3*(nParents-1)*nSets);              	% Upper bounds of parameters
            % Randomly set the parameters
            obj.parameter = unifrnd(obj.lower,obj.upper);
            obj.ParameterAssign();
        end
        %% Assign parameters to variables
        function ParameterAssign(obj)
            obj.Weight = reshape(obj.parameter,[],obj.nSets)';
            obj.Fit    = cumsum(obj.Weight(:,3:3:end),1);
            obj.Fit    = obj.Fit./repmat(max(obj.Fit,[],1),size(obj.Fit,1),1);
        end
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            ParentDec = obj.Gather(Problem,Precursors,Ratio,2,obj.nParents);
            R = zeros(size(ParentDec,1)-size(ParentDec,1)./obj.nParents,size(ParentDec,2));
            for i = 1 : obj.nParents-1
                R(i:obj.nParents-1:end,:) = ParaSampling([size(ParentDec,1)./obj.nParents,size(R,2)],obj.Weight(:,(i-1)*3+1:(i-1)*3+2),obj.Fit(:,i));
            end
            OffDec = zeros(size(ParentDec,1)./obj.nParents,size(ParentDec,2));
            for i = 1 : size(OffDec,1)
                r = R((i-1)*(obj.nParents-1)+1:i*(obj.nParents-1),:);
                OffDec(i,:) = sum([1-sum(r,1);r].*ParentDec((i-1)*obj.nParents+1:i*obj.nParents,:),1);
            end
            obj.output = min(repmat(Problem.upper,size(OffDec,1),1),max(repmat(Problem.lower,size(OffDec,1),1),OffDec));
        end
    end
end

function r = ParaSampling(xy,weight,fit)
% Parameter sampling

    r    = repmat(randn(xy(1),1),1,xy(2));
    type = arrayfun(@(S)find(rand<=fit,1),zeros(xy));
    for i = 1 : length(fit)
        index = type == i;
        r(index) = r(index)*weight(i,1) + weight(i,2);
    end
end