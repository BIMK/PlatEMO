classdef Block_Kopt < BLOCK
% k-opt
% k --- 4 --- Max number of k for k-opt

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        k;      % <hyperparameter> Max number of k for k-opt
        Fit;	% <parameter> Probability of using k-opt
    end
    methods
        %% Default settings of the block
        function obj = Block_Kopt(k)
            obj.k     = k;
            obj.lower = zeros(1,obj.k) + 1e-20;
            obj.upper = ones(1,obj.k);
            % Randomly set the parameters
            obj.parameter = unifrnd(obj.lower,obj.upper);
            obj.ParameterAssign();
        end
        %% Assign parameters to variables
        function ParameterAssign(obj,~,~)
            obj.Fit = cumsum(obj.parameter);
            obj.Fit = obj.Fit./max(obj.Fit);
        end
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            ParentDec = obj.Gather(Problem,Precursors,Ratio,2,1);
            type = arrayfun(@(S)find(rand<=obj.Fit,1),1:size(ParentDec,1));
            for i = find(type>1)
                s = randi([1 Problem.D-(type(i)-1)*2],1,type(i));
                for j = 2 : type(i)
                    s(j) = randi([s(j-1)+2,Problem.D-(type(i)-j)*2]);
                end
                newPerm = 1 : s(1);
                for j = randperm(type(i)-1)
                    if type(i)==2 || rand>0.5
                        newPerm = [newPerm,flip(s(j)+1:s(j+1))];
                    else
                        newPerm = [newPerm,s(j)+1:s(j+1)];
                    end
                end
                newPerm = [newPerm,s(end)+1:Problem.D];
                ParentDec(i,:) = ParentDec(i,newPerm);
            end
            obj.output = ParentDec;
        end
    end
end