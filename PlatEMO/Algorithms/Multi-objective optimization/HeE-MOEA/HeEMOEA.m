classdef HeEMOEA < ALGORITHM
% <multi> <real/integer> <expensive>
% Multiobjective evolutionary algorithm with heterogeneous ensemble based
% infill criterion
% Ke --- 5 --- Number of the solutions to be revaluated

%------------------------------- Reference --------------------------------
% D. Guo, Y. Jin, J. Ding, and T. Chai. Heterogeneous ensemble-based infill 
% criterion for evolutionary multiobjective optimization of expensive 
% problems. IEEE Transactions on Cybernetics, 2019, 49(3): 1012-1025.
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
            assert(~isempty(ver('nnet')),'The execution of HeE-MOEA requires the Deep Learning Toolbox.');
            
            %% Parameter setting
            Ke = Algorithm.ParameterSet(5);
            NI = 11*Problem.D-1;
            P  = UniformPoint(NI,Problem.D,'Latin');
            A  = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            
            %% Settings of the Ensemble Model
            L    = 11*Problem.D-1+25;
            ADec = A.decs;
            AObj = A.objs;
            
            Selection = PSOMyself(ADec, AObj(:,1),Problem.D);
            str1={'FE', 'FS', 'NONE'};
            str2={'RBF1','SVM', 'RBF2'};%, 'KNN', 'DTree'
            for i=1:length(str1)
                for j=1:length(str2)
                    str{(i-1)*length(str2)+j, 1}=str1{i};
                    str{(i-1)*length(str2)+j, 2}=str2{j};
                end
            end

            while Algorithm.NotTerminated(A)                
                %% Update the model 
              % Select the train data
                ADec = A.decs;
                AObj = A.objs;
                Numdata = size(ADec,1);
                if Numdata <=L
                    % fprintf('No training data decrease\n');                    
                else
                    FrontNo   = NDSort(AObj,Numdata);
                    [~,index] = sort(FrontNo);
                    ADec1 = ADec(index(1:floor(L/2)), :);
                    AObj1 = AObj(index(1:floor(L/2)), :);
                    ADec2 = ADec(index(floor(L/2)+1:end), :);
                    AObj2 = AObj(index(floor(L/2)+1:end), :);
                    index = randperm(size(ADec2,1));
                    ADec  = [ADec1;ADec2(index(1:L-floor(L/2)),:)];
                    AObj  = [AObj1;AObj2(index(1:L-floor(L/2)),:)];
                end
              
                % Train the model
                Models = TrainModel(ADec,AObj,Selection,str,Problem.M,Problem.D);
              
                % Optimization
                New    = NSGA2ESelection(ADec,AObj,Models,str,Problem,Ke);
                PopNew = Problem.Evaluation(New);
                A      = [A PopNew];
            end
        end
    end
end