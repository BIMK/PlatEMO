classdef SLMEA < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Super-large-scale multi-objective evolutionary algorithm
% useGPU --- 0 --- Whether use GPU acceleration

%------------------------------- Reference --------------------------------
% Y. Tian, Y. Feng, X. Zhang, and C. Sun, A fast clustering based
% evolutionary algorithm for super-large-scale sparse multi-objective
% optimization, IEEE/CAA Journal of Automatica Sinica, 2021.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

     methods
         function main(Algorithm,Problem)
            %% Parameter setting
            useGPU = Algorithm.ParameterSet(0);
             
            %% Generate initial variables
            REAL        = ~strcmp(Problem.encoding,'binary');
            P           = 0.5;
            ArchiveObj  = [];
            ArchiveMask = [];
            numGroup    = 10;
            Lr          = 0;
            RC          = 1;
            T           = 20;
            GlobalGen   = 1;
            
            %% Population initialization
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),ones(1,Problem.D))) = 1;
            end
            if ~useGPU
                Lower = Problem.lower;
                Upper = Problem.upper;
                if REAL
                    Dec = unifrnd(repmat(Lower,Problem.N,1),repmat(Upper,Problem.N,1));
                else
                    Dec = ones(Problem.N,Problem.D);
                end
            else
                Lower = gpuArray(Problem.lower);
                Upper = gpuArray(Problem.upper);
                Mask  = gpuArray(Mask);
                if REAL
                    Dec = unifrnd(repmat(Lower,Problem.N,1),repmat(Upper,Problem.N,1));
                else
                    Dec = gpuArray.ones(Problem.N,Problem.D);
                end
            end            
            Population = SOLUTION(Dec.*Mask);
            
           %% Update Archive
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N);
            PopObj      = gather(Population.objs);
            ArchiveObj   = [ArchiveObj;PopObj(FrontNo==1,:)];
            ArchiveMask = [ArchiveMask;Mask(FrontNo==1,:)];
            [Fitness,Mix,Zero,One] = CalculateFitness(ArchiveMask,Problem.D);
            
           %% optimization
            while Algorithm.NotTerminated(Population)
                % Update population
                MatingPool = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask,long,numGroup] = Operator(Dec(MatingPool,:),Mask(MatingPool,:),REAL,Fitness,Mix,P,T,Zero,One,Upper,Lower,RC,numGroup,useGPU,GlobalGen);
                Offspring  = SOLUTION(OffDec.*OffMask);
                GlobalGen  = GlobalGen + 1;
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                % Update P
                if GlobalGen >= T+1
                  Fn  = NDSort(gather(Offspring.objs),1); 
                  LOC = find(Fn==1);
                  a   = find(LOC<=long);
                  b   = find(LOC>long);
                  P   = min(0.95,max(0.05,0.5*(P+(length(a)*(Problem.N-long)/(length(a)*(Problem.N-long)+length(b)*long)))));
                  RC  = exp((length(a)/(long+0.00001)-Lr)/numGroup);
                  Lr  = length(a)/long;
                end
                % Update Archive
                PopObj          = gather(Population.objs);
                ArchiveObj      = [ArchiveObj;PopObj(FrontNo==1,:)];
                ArchiveMask     = [ArchiveMask;Mask(FrontNo==1,:)];
                [ArchiveObj,ia] = unique(ArchiveObj,'rows');
                [frontNo,~]     = NDSort(ArchiveObj,1);
                ArchiveObj      = ArchiveObj(frontNo==1,:);
                ArchiveMask     = ArchiveMask(ia,:);
                ArchiveMask     = ArchiveMask(frontNo==1,:);
                if size(ArchiveObj,1) > 100
                    ArLocation  = randperm(size(ArchiveObj,1),100);
                    ArchiveObj  = ArchiveObj(ArLocation,:);
                    ArchiveMask = ArchiveMask(ArLocation,:);
                end
                % Calculate Fitness
                [Fitness,Mix,Zero,One] = CalculateFitness(ArchiveMask,Problem.D);
             end
         end
     end
end