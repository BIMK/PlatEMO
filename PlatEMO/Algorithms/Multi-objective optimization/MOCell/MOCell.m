classdef MOCell < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% Cellular genetic algorithm

%------------------------------- Reference --------------------------------
% A. J. Nebro, J. J. Durillo, F. Luna, B. Dorronsoro, and E. Alba, MOCell:
% A cellular genetic algorithm for multiobjective optimization,
% International Journal of Intelligent Systems, 2009, 24(7): 726-746.
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
            %% Neighborhood of each cell
            cSize    = floor(sqrt(Problem.N));
            n1       = reshape(1:cSize^2,cSize,cSize);
            n2       = n1([end,1:end-1],[end,1:end-1]);
            n3       = n1([end,1:end-1],:);
            n4       = n1([end,1:end-1],[2:end,1]);
            n5       = n1(:,[end,1:end-1]);
            n6       = n1(:,[2:end,1]);
            n7       = n1([2:end,1],[end,1:end-1]);
            n8       = n1([2:end,1],:);
            n9       = n1([2:end,1],[2:end,1]);
            Neighbor = [n1(:),n2(:),n3(:),n4(:),n5(:),n6(:),n7(:),n8(:),n9(:)];

            %% Generate random population
            Population = Problem.Initialization(size(Neighbor,1));
            Archive    = [];

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                % Update the fitness of solutions in the population
                FrontNo  = NDSort(Population.objs,Population.cons,inf);
                CrowdDis = CrowdingDistance(Population.objs,FrontNo);

                % Generate one offspring for each cell
                Offspring(1:length(Population)) = SOLUTION();
                for i = 1 : length(Population)
                    parents      = Neighbor(i,TournamentSelection(2,2,FrontNo(Neighbor(i,:)),-CrowdDis(Neighbor(i,:))));
                    Offspring(i) = OperatorGAhalf(Problem,Population(parents));
                end

                % Replace the old population
                ViolationP = sum(max(0,Population.cons),2);
                ViolationO = sum(max(0,Offspring.cons),2);
                Replace    = ViolationO < ViolationP | ViolationO == ViolationP & any(Population.objs>Offspring.objs,2);
                Population(Replace) = Offspring(Replace);

                % Update the archive
                Archive = [Archive,Offspring];
                Archive = Archive(NDSort(Archive.objs,Archive.cons,1)==1);
                if length(Archive) > Problem.N
                    [~,rank] = sort(CrowdingDistance(Archive.objs,ones(1,length(Archive))),'descend');
                    Archive  = Archive(rank(1:Problem.N));
                end

                % Feedback
                nRep = min(min(20,length(Population)),length(Archive));
                Population(randperm(length(Population),nRep)) = Archive(randperm(length(Archive),nRep));
            end
        end
    end
end