classdef EDNARMOEA < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Efficient dropout neural network based AR-MOEA
% delta --- 0.05 --- Threshold of judging the diversity
% wmax  ---   20 --- Number of generations before updating Kriging models
% Ke    ---    3 --- Number of the solutions to be revaluated in each iteration

%------------------------------- Reference --------------------------------
% D. Guo, X. Wang, K. Gao, Y. Jin, J. Ding, and T. Chai. Evolutionary
% optimization of high-dimensional multiobjective and many-objective
% expensive problems assisted by a dropout neural network. IEEE
% Transactions on Systems, Man, and Cybernetics: Systems, 2022, 52(4):
% 2084-2097.
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
            assert(~isempty(ver('nnet')),'The execution of EDN-ARMOEA requires the Deep Learning Toolbox.');
            
            %% Parameter setting
            [delta,wmax,Ke] = Algorithm.ParameterSet(0.05,20,3);            
            W     = UniformPoint(Problem.N,Problem.M);
            NI    = 11*Problem.D-1;
            P     = UniformPoint(NI,Problem.D,'Latin');
            A     = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            tr_x  = A.decs;
            tr_y  = A.objs;
            [tr_xx,ps] = mapminmax(tr_x');tr_xx=tr_xx';
            [tr_yy,qs] = mapminmax(tr_y');tr_yy=tr_yy';
            Params.ps  = ps;Params.qs=qs;           
            RatioOld   = [];
            
            %% Train the model
            [net, Params] = trainmodel(tr_xx, tr_yy, Params);

            while Algorithm.NotTerminated(A)
                %% Update the model 
                net=updatemodel(tr_xx, tr_yy, Params, net);

                %% Generate the sampling points and random population
                popsize=Problem.N;
                PopDec=repmat(Problem.upper-Problem.lower,popsize,1).*rand(popsize,Problem.D)+repmat(Problem.lower,popsize,1);
                [PopObj, PopMSE]=Estimate(PopDec, net, Params, Problem.M);
                [Archive,RefPoint,Range, Ratio] = UpdateRefPoint(PopObj,W,[]);
                if isempty(RatioOld)
                    RatioOld=Ratio;
                end
               
                %% Start the interations
                w=1;
                while w < wmax
                    MatingPool = MatingSelection(PopObj,RefPoint,Range);
                    OffspringDec  = OperatorGA(Problem,PopDec(MatingPool,:),{1,20,1,20});
                    [OffspringObj, OffspringMSE] = Estimate(OffspringDec, net, Params, Problem.M);
                    [Archive,RefPoint,Range, Ratio] = UpdateRefPoint([Archive;OffspringObj],W,Range);
                    MediatePopDec=[PopDec;OffspringDec];
                    MediatePopObj=[PopObj;OffspringObj];
                    MediatePopMSE=[PopMSE;OffspringMSE];
                    [Index,Range]       = EnvironmentalSelection(MediatePopObj,RefPoint,Range,popsize);
                    PopDec=MediatePopDec(Index,:);
                    PopObj=MediatePopObj(Index,:);
                    PopMSE=MediatePopMSE(Index,:);
                    w=w+1;
                end
                flag=RatioOld-Ratio<delta;
                PopNew=IndividualSelect(PopDec, PopObj, PopMSE, Ke, flag);
                RatioOld=Ratio;
                New = Problem.Evaluation(PopNew);
                A   = [A,New];
                [tr_x, tr_y]=SelectTrainData(A, 11*Problem.D-1, length(New));
                [tr_xx,ps]=mapminmax(tr_x');tr_xx=tr_xx';
                [tr_yy,qs]=mapminmax(tr_y');tr_yy=tr_yy';
                Params.ps=ps;Params.qs=qs;   
            end
        end
    end
end