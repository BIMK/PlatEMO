classdef Block_Mutation < BLOCK
% Unified mutation for real variables
% nSets --- 5 --- Number of parameter sets

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nSets;      % <hyperparameter> Number of weight sets
        Weight;     % <parameter> Weight sets
        Fit;        % <parameter> Expectation of using each weight set
        nDec = 1;   % Number of decision variables
    end
    methods
        %% Default settings of the block
        function obj = Block_Mutation(nSets)
            obj.nSets = nSets;	% Hyperparameter
            obj.lower = repmat([0 1e-20],1,nSets);	% Lower bounds of parameters
            obj.upper = repmat([1 5],1,nSets);  	% Upper bounds of parameters
            % Randomly set the parameters
            obj.parameter = unifrnd(obj.lower,ones(1,2*nSets));
            obj.ParameterAssign();
        end
        %% Assign parameters to variables
        function ParameterAssign(obj)
            obj.Weight = reshape(obj.parameter,[],obj.nSets)';
            obj.Weight(:,end) = obj.Weight(:,end)./obj.nDec;
            obj.Weight = [obj.Weight;0,max(0,1-sum(obj.Weight(:,end)))];
            obj.Fit    = cumsum(obj.Weight(:,end));
            obj.Fit    = obj.Fit./max(obj.Fit);
        end
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            ParentDec = obj.Gather(Problem,Precursors,Ratio,2,1);
            if size(ParentDec,2) ~= obj.nDec
                obj.nDec = size(ParentDec,2);
                obj.ParameterAssign();
            end
            r          = ParaSampling(size(ParentDec),obj.Weight(:,1),obj.Fit);
            obj.output = ParentDec + repmat(Problem.upper-Problem.lower,size(ParentDec,1),1).*r;
        end
    end
end

function r = ParaSampling(xy,weight,fit)
% Parameter sampling

    r    = repmat(randn(xy(1),1),1,xy(2));
    type = arrayfun(@(S)find(rand<=fit,1),zeros(xy));
    for i = 1 : length(fit)
        index = type == i;
        r(index) = r(index)*weight(i,1);
    end
end