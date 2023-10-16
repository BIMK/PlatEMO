classdef REMO < ALGORITHM
% <multi/many> <real> <expensive>
% Expensive multiobjective optimization by relation learning and prediction
% k    ---    6 --- Number of reference solutions
% gmax --- 3000 --- Number of solutions evaluated by surrogate model

%------------------------------- Reference --------------------------------
% H. Hao, A. Zhou, H. Qian, and H. Zhang, Expensive multiobjective
% optimization by relation learning and prediction, IEEE Transactions on
% Evolutionary Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameterr setting
            [k,gmax] = Algorithm.ParameterSet(6,3000);

            %% Initalize the population by Latin hypercube sampling
            if Problem.D <= 10
                N = 11*Problem.D-1;
            else
                N = 100;
            end
            PopDec     = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            Archive    = Population;

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                % Select reference solutions and preprocess the data
                Ref       = RefSelect(Population,k);
                Input     = Population.decs; 
                Catalog   = GetOutput_PBI(Population.objs,Ref.objs); 
                [XXs,YYs] = GetRelationPairs(Input,Catalog);
                [TrainIn,TrainOut,TestIn,TestOut] = DataProcess(XXs,YYs);
                xDim = size(TrainIn,2);
                
                % Train relation model
                [TrainIn_nor,TrainIn_struct] = mapminmax(TrainIn');
                TrainIn_nor     = TrainIn_nor';
                TrainOut_onehot = onehotconv(TrainOut,1);
                net = patternnet([ceil(xDim*1.5),xDim*1,ceil(xDim/2)]);
                net.trainParam.showWindow =0;
                net        = train(net,TrainIn_nor',TrainOut_onehot');
                TestIn_nor = mapminmax('apply',TestIn',TrainIn_struct)';
                TestPre    = onehotconv(net(TestIn_nor')',2);             
                p_err      = sum(TestPre ~= TestOut)/size(TestPre,1);
                Smodel.X   = Input;
                Smodel.Y   = Catalog;
                Smodel.mp_struct = TrainIn_struct;
                Smodel.net       = net;
                Smodel.p_err     = p_err;
                Next = RSurrogateAssistedSelection(Problem,Ref,Population.decs,gmax,Smodel);
                if ~isempty(Next)
                    Archive = [Archive,Problem.Evaluation(Next)];
                end
                Population = RefSelect(Archive,Problem.N);
            end
        end
	end
end