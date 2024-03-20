classdef MOEADD < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Many-objective evolutionary algorithm based on dominance and
% decomposition
% delta --- 0.9 --- The probability of choosing parents locally

%------------------------------- Reference --------------------------------
% K. Li, K. Deb, Q. Zhang, and S. Kwong, An evolutionary many-objective
% optimization algorithm based on dominance and decomposition, IEEE
% Transactions Evolutionary Computation, 2015, 19(5): 694-716.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            delta = Algorithm.ParameterSet(0.9);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            [~,Region] = max(1-pdist2(Population.objs,W,'cosine'),[],2);
            FrontNo    = NDSort(Population.objs,inf);
            Z          = min(Population.objs,[],1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                for i = 1 : Problem.N            
                    % Choose the parents
                    Ei = find(ismember(Region,B(i,:)));
                    if rand < delta && length(Ei) >= 2
                        P = Ei(TournamentSelection(2,2,sum(max(0,Population(Ei).cons),2)));
                    else
                        P = TournamentSelection(2,2,sum(max(0,Population.cons),2));
                    end

                    % Generate an offspring
                    Offspring = OperatorGAhalf(Problem,Population(P));
                    [~,offRegion] = max(1-pdist2(Offspring.obj,W,'cosine'));

                    % Add the offspring to the population
                    Population = [Population,Offspring];
                    PopObj     = Population.objs;
                    Region     = [Region;offRegion];
                    FrontNo    = UpdateFront(PopObj,FrontNo);
                    CV         = sum(max(0,Population.cons),2);

                    % Update the ideal point
                    Z = min(Z,Offspring.obj);

                    % Delete a solution from the population
                    if any(CV>0)
                        [~,S] = sort(CV,'descend');
                        S     = S(1:sum(CV>0));
                        flag  = false;
                        for j = 1 : length(S)
                            if sum(Region==Region(S(j))) > 1
                                flag = true;
                                x    = S(j);
                                break;
                            end
                        end
                        if ~flag
                            x = S(1);
                        end
                    elseif max(FrontNo) == 1
                        x = LocateWorst(PopObj,W,Region,FrontNo,Z);
                    else
                        Fl = find(FrontNo==max(FrontNo));
                        if length(Fl) == 1
                            if sum(Region==Region(Fl)) > 1
                                x = Fl;
                            else
                                x = LocateWorst(PopObj,W,Region,FrontNo,Z);
                            end
                        else
                            SubRegion = unique(Region(Fl));
                            Crowd     = hist(Region(ismember(Region,SubRegion)),1:size(W,1));
                            Phi       = find(Crowd==max(Crowd));
                            PBI       = CalPBI(PopObj,W,Region,Z,ismember(Region,Phi));
                            PBISum    = zeros(1,size(W,1));
                            for j = 1 : length(PBI)
                                PBISum(Region(j)) = PBISum(Region(j)) + PBI(j);
                            end
                            [~,Phi] = max(PBISum);
                            Phih    = find(Region==Phi);
                            if length(Phih) > 1
                                [~,x] = max(PBI(Phih));
                                x = Phih(x);
                            else
                                x = LocateWorst(PopObj,W,Region,FrontNo,Z);
                            end
                        end
                    end
                    Population(x) = [];
                    Region(x)     = [];
                    FrontNo       = UpdateFront(Population.objs,FrontNo,x);
                end
            end
        end
    end
end