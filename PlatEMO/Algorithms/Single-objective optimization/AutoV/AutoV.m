classdef AutoV < ALGORITHM
% <2023> <single/multi> <real/integer> <large/none> <constrained/none>
% Automated design of variation operators
% fileName --- 'Weight.mat' --- File name for saving/loading weights
% Mode     ---            1 --- 1. Test 2. Training
% TrainN   ---           50 --- Population size of GA for training
% TrainFE  ---         5000 --- Maximum evaluations of GA for training

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, C. He, K. C. Tan, and Y. Jin. Principled design of
% translation, scale, and rotation invariant variation operators for
% metaheuristics. Chinese Journal of Electronics, 2023, 32(1): 111-129.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            [fileName,Mode,TrainN,TrainFE] = Algorithm.ParameterSet('Weight.mat',1,50,5000);
            if Mode == 1    % Test mode
                % Load weights
                if ischar(fileName)
                    try
                        load(fullfile(fileparts(mfilename('fullpath')),fileName),'-mat','Population');
                    catch
                        error('Fail to load weights from %s. This algorithm should be trained before used.',fileName);
                    end
                    Weight = reshape(Population(1).dec,[],4);
                else
                    Weight = fileName;
                end
                Fit = cumsum(Weight(:,end));
                Fit = Fit./max(Fit);
                % Optimization
                Population = Problem.Initialization();
                [Population,Fitness] = EnvironmentalSelection(Population,Problem.N);
                while Algorithm.NotTerminated(Population)
                    MatingPool = TournamentSelection(2,2*Problem.N,Fitness);
                    Offspring  = TSRIOperator(Problem,Weight,Fit,Population(MatingPool));
                    [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N);
                end
            else            % Training mode
                Lower = reshape(repmat([0 0 -1 1e-6],10,1),1,[]);
                Upper = reshape(repmat([1 1  1    1],10,1),1,[]);
                platemo('algorithm',@GA,'outputFcn',@trainOutputFcn,'initFcn',@trainInit,'objFcn',@trainObj,...
                        'data',{Problem,fullfile(fileparts(mfilename('fullpath')),fileName),Lower,Upper},...
                        'once',true,'N',TrainN,'maxFE',TrainFE,'D',length(Lower),'lower',Lower,'upper',Upper);
            end
        end
    end
end

function PopObj = trainObj(PopDec,data)
% Objective function for training

    PopObj = zeros(size(PopDec,1),3);
    for i = 1 : numel(PopObj)
        Future(i) = parfeval(@parallelFcn,1,data{1},PopDec(mod(i-1,size(PopDec,1))+1,:));
    end
    while ~all([Future.Read])
        [r,Population] = fetchNext(Future,0.01);
        if ~isempty(r)
            [i,j] = ind2sub(size(PopObj),r);
            if data{1}.M == 1
                PopObj(i,j) = data{1}.CalMetric('Min_value',Population);
            else
                PopObj(i,j) = -data{1}.CalMetric('HV',Population);
            end
        end
    end
    PopObj = median(PopObj,2);  % Median best value among 3 runs
end

function Population = parallelFcn(Problem,Dec)
% Evaluation of a single operator

    Algorithm = AutoV('parameter',{reshape(Dec,[],4)},'outputFcn',@(~,~)[]);
    Algorithm.Solve(Problem);
    Population = Algorithm.result{end};
end

function PopDec = trainInit(N,data)
% Initialization function for training

    if exist(data{2},'file')
        load(data{2},'-mat','Population');
        PopDec = Population.decs;
    else
        PopDec = [];
    end
    if size(PopDec,1) > N
        PopDec = PopDec(1:N,:);
    elseif size(PopDec,1) < N
        PopDec = [PopDec;unifrnd(repmat(data{3},N-size(PopDec,1),1),repmat(data{4},N-size(PopDec,1),1))];
    end
end

function trainOutputFcn(ALG,PRO)
% Output function for training

    Population = ALG.result{end};
    clc; fprintf('Training on %s for %d/%d evaluations, %.1fs passed, the best value is %.4e\n',class(PRO.data{1}),PRO.FE,PRO.maxFE,ALG.metric.runtime,min(Population.objs));
    save(PRO.data{2},'Population');
end