classdef MOEACKF < ALGORITHM
% <2024> <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Multi-objective evolutionary algorithm based on cross-scale knowledge fusion

%------------------------------- Reference --------------------------------
% Z. Ding, L. Chen, D. Sun, and X. Zhang. Efficient sparse large-scale
% multi-objective optimization based on cross-scale knowledge fusion. IEEE
% Transactions on Systems, Man, and Cybernetics: Systems, 2024, 54(11):
% 6989-7001.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lei Chen (email: soarfx@163.com)

    methods
        function main(Algorithm,Problem)
            %% Population initialization
            REAL = ~strcmp(Problem.encoding,'binary');
            [Dec,Mask,Population,Fitness]=PriorAnalysis_initialization(Problem,REAL);
            [Population,Dec,Mask] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);
            rho   = 0.5;
            GROUP = [];     % In order to prevent outliers in Kmeans algorithm
           
            %% Optimization
            while Algorithm.NotTerminated(Population)  
                [FrontNo,~] = NDSort(Population.objs,Population.cons,inf);
                CrowdDis    = CrowdingDistance(Population.objs,FrontNo);
                % Divide the population, Pop1 is used as the population for 
                % dual dimension reduction to generate offspring, and Pop2 
                % is used as the population of dual grouping.
                Pop1_Site = rho > rand(1,length(Population));                
                if sum(Pop1_Site) >= 2
                    Pop1  = Population(Pop1_Site);
                    Dec1  = Dec(Pop1_Site,:);
                    Mask1 = Mask(Pop1_Site,:);
                    Pop2  = Population(~Pop1_Site);
                    Dec2  = Dec(~Pop1_Site,:);
                    Mask2 = Mask(~Pop1_Site,:);
                elseif sum(~Pop1_Site) < 1
                    Pop1  = Population;
                    Dec1  = Dec;
                    Mask1 = Mask;
                    [Pop2,Dec2,Mask2] = deal([]);                
                else
                    [Pop1,Dec1,Mask1] = deal([]);
                    Pop2  = Population;
                    Dec2  = Dec;
                    Mask2 = Mask;     
                end
                % decision variable Sparsity analysis
                [Local_Knowlege,Global_Knowledge,Fitness,NSV,SV,theta] = SparsityAnalysis(Problem,Mask,FrontNo,Fitness,GROUP);
                GROUP = [NSV;SV];
                % dual dimension reduction
               if ~isempty(Pop1)
                   [OffDec1,OffMask1,len1] = Reproduction1(Problem,Pop1,Dec1,Mask1,FrontNo,CrowdDis,Pop1_Site,Local_Knowlege,Global_Knowledge,NSV,SV,theta,REAL);
                   Offspring1 = Problem.Evaluation(OffDec1.*OffMask1);
               else
                   [OffDec1,OffMask1,Offspring1] = deal([]);
                   len1 = 0;
               end
               % dual grouping 
               if ~isempty(Pop2)
                   MatingPool = TournamentSelection(2,(Problem.N-len1)*2,FrontNo(~Pop1_Site),-CrowdDis(~Pop1_Site));
                   [OffDec2,OffMask2] = Reproduction2(Problem,Dec2(MatingPool,:),Mask2(MatingPool,:),SV,NSV,REAL);
                   Offspring2 = Problem.Evaluation(OffDec2.*OffMask2);
               else
                  [OffDec2,OffMask2,Offspring2] = deal([]);
               end
               [Population,Dec,Mask,sRatio] = EnvironmentalSelection([Population,Offspring1,Offspring2],[Dec;OffDec1;OffDec2],[Mask;OffMask1;OffMask2],Problem.N,length(Population),len1);
               rho = (rho+sRatio)/2;
            end
        end
    end
end