classdef AdaW < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Evolutionary algorithm with adaptive weights

%------------------------------- Reference --------------------------------
% M. Li and X. Yao, What weights work for you? Adapting weights for any
% Pareto front shape in decomposition-based evolutionary multiobjective
% optimisation, Evolutionary Computation, 2020, 28(2): 227-253.
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
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);  % Generate the weight vectors
            T = ceil(Problem.N/10);                             % The size of neighbours of each weight
            
            %% Detect the neighbours of each weight
            B = pdist2(W,W);
            [~,B] = sort(B,2); 
            B = B(:,1:T);
            
            %% Generate random population
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);
            
            %% Generate an archive set
            Archive = Population(NDSort(Population.objs,1)==1);
            Archive_temp = Population; 
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                  % For each weight
                  for i = 1 : Problem.N 
                      % Choose parents
                      if rand < 0.9
                          P = B(i,randperm(size(B,2)));
                      else
                          P = randperm(Problem.N);
                      end
                      % Generate an offspring
                      Offspring = OperatorGAhalf(Problem,Population(P(1:2)));
                      % Put the offspring into the archive
                      Archive_temp(i) = Offspring;
                      % Update the ideal point
                      Z = min(Z ,Offspring.obj);
                      % Pick a neighbour to update   
                      g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
                      g_new = max(repmat(abs(Offspring.obj-Z),length(P),1)./W(P,:),[],2);
                      Population(P(find(g_old >= g_new,1))) = Offspring;
                  end
                  Archive = [Archive,Archive_temp];
                  % Maintenance operation in the archive set
                  Archive = Archive(NDSort(Archive.objs,1)==1);
                  Archive = ArchiveUpdate(Archive, 2 * Problem.N);
                  % Update weights
                  if ~mod(ceil(Problem.FE/Problem.N),ceil(0.05*ceil(Problem.maxFE/Problem.N))) && Problem.FE <= Problem.maxFE*0.9
                     [Population,W,B] = WeightUpdate(Population,W,Archive,Z,T,Problem);
                  end
            end
        end
    end
end