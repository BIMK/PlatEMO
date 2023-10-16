classdef Sparse_PM < PROBLEM
% <multi> <binary> <large/none> <expensive/none> <sparse/none>
% The pattern mining problem
% numTran --- 10000 --- Number of transactions
% lenTran ---    50 --- Average length of transactions
% numPa   ---   100 --- Number of patterns
% lenPa   ---     5 --- Average length of patterns
% numItem ---   100 --- Number of items (decision variables)

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, C. Wang, and Y. Jin, An evolutionary algorithm for
% large-scale sparse multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2020, 24(2): 380-393.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Data;   % Transaction dataset
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [numTran,lenTran,numPa,lenPa,numItem] = obj.ParameterSet(10000,50,100,5,100);
            obj.M        = 2;
            obj.D        = numItem;
            obj.encoding = 4 + zeros(1,obj.D);
            % Randomly generate transaction dataset
            fileName = fullfile(fileparts(mfilename('fullpath')),sprintf('Dataset_PM-D%d-T%d-L%d-I%d-N%d.mat',numTran,lenTran,numPa,lenPa,numItem));
            if exist(fileName,'file') == 2
                load(fileName,'Data');
            else
                % Generate patterns
                PaSet = cell(1,numPa);
                RAND1 = min(numItem,max(1,random('Poisson',lenPa,1,numPa)));    % Length of each pattern
                RAND2 = min(1,max(0,random('Exponential',0.5,1,numPa)));	% Fraction of items in each pattern chosen from previous pattern
                PaSet{1} = randperm(numItem,RAND1(1));
                for i = 2 : numPa
                    PaSet{i} = [PaSet{i-1}(randperm(end,min(end,round(RAND2(i)*RAND1(i))))),...
                                randperm(numItem,round((1-RAND2(i))*RAND1(i)))];
                end
                % Generate transactions
                RAND3 = max(0,random('Exponential',1,1,numPa));	% Probability that each pattern will be picked
                RAND3 = cumsum(RAND3)./sum(RAND3);
                RAND4 = random('Normal',0.5,0.1,1,numPa);	% Corruption level of each pattern
                Data  = false(numTran,numItem);
                RAND5 = min(min(numItem,length(unique(cat(2,PaSet{:})))),max(1,random('Poisson',lenTran,1,numTran)));	% Length of each transaction
                Data = false(numTran,numItem);
                for i = 1 : numTran
                    while sum(Data(i,:)) < RAND5(i)
                        current = find(rand<=RAND3,1);
                        Data(i,PaSet{current}(rand(1,end)>RAND4(current))) = true;
                    end
                end
                save(fileName,'Data');
            end
            obj.Data = Data;
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2
                N = obj.N;
            end
            PopDec     = rand(N,obj.D) < 20/obj.D;
            Population = obj.Evaluation(PopDec);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(PopDec);
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                x  = PopDec(i,:);
                Tx = all(obj.Data(:,x),2);
                if ~any(Tx)
                    PopObj(i,:) = 1;
                else
                    PopObj(i,1) = 1 - mean(Tx);
                    PopObj(i,2) = 1 - mean(sum(x)./sum(obj.Data(Tx,:),2));
                end
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(1-Population.objs,{'Frequency','Occupancy',[]});
        end
    end
end