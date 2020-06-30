function S3CMAES(Global)
% <algorithm> <S>
% Scalable small subpopulations based covariance matrix adaptation
% evolution strategy

%------------------------------- Reference --------------------------------
% H. Chen, R. Cheng, J. Wen, H. Li, and J. Weng, Solving large-scale
% many-objective optimization problems by covariance matrix adaptation
% evolution strategy with scalable small subpopulations, Information
% Sciences, 2018.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen
    
    %% Detect the group of each distance variable
    nPer      = 50;	% Sample size to divide the convergence- and diversity-related variables
    nPerGroup = 35;	% the group size for separative variables
    
    [PV, DV] = ControlVariableAnalysis(Global,nPer);	% divide the convergence- and diversity-related variables        
    Groups   = GroupDV(Global,DV,PV,nPerGroup);         % divide the convergence-related variables based on correlation
    
    popN = 5;	% the number of sub-populations
    % V: random unit vectors in diversity space
    V = 0.05 + 0.9*rand(popN,length(PV));
    % extend to the diversity space
    V = repmat(Global.lower(PV),popN,1)+V.*(repmat((Global.upper(PV)-Global.lower(PV)),popN,1));  
    
    % CMA-ES parameters
    popSize         = 6 + floor(3*log(nPerGroup));	% the population size for CMA-ES
    tempParaOnePopu = LoadCMAESparameters(Global,Groups,popSize);
    CMAParaMPopu    = repmat(tempParaOnePopu,popN,1);
    
    BigPopulation = cell(1,popN);	% initialize the big population to contain all the population
    % the diversity-related variables of all the solutions in a sub-population are the same
    for i = 1: popN
        tempDecs = zeros(popSize,Global.D);
        tempDecs(:, PV)  = repmat(V(i, :),popSize,1);
        DVPositionM      = Global.lower(DV) + (Global.upper(DV)-Global.lower(DV)).*rand(1, length(DV));
        tempDecs(:, DV)  = repmat(DVPositionM,popSize,1);
        BigPopulation{i} = tempDecs;
    end
    
    %% Optimization
    Archive      = INDIVIDUAL(BigPopulation{1});
    lastBestValA = zeros(1,popN);	% Record the best val last iteration
    stopTag      = false(1,popN);
    unUpdateNum  = zeros(1,popN);
    unChangeThr  = 1e-10;
    firstTag     = true;
    
    ConvergedSolutionSet = [];	% Record the converged solutions
    while Global.NotTermination(Archive)
        
        Archive = ConvergedSolutionSet;
        popN    = length(BigPopulation);
        for p = 1: popN	% evolute each small population
            if stopTag(p)	% if the subpopulaton has converged, then no evolve it
                continue;
            end
            
            tempDecs = BigPopulation{p};	% obtain the p-th population
            PVDecs   = tempDecs(:,PV);
            bestmem  = tempDecs(1,:);   	% select the first individual as best member
            
            tempDecs        = zeros(popSize,Global.D);	% record the new population
            tempDecs(:, PV) = PVDecs;                   % repmat(V(p, :), popSize, 1);
            
            for g = 1:  length(Groups) % evolute each group of convergence-related variables for a sub-population
                dim_index = Groups{g};
                % employ the CMA-ES
                [CMAParaMPopu(p,g),pop,BestVal,BestIndividual] = CMAES(CMAParaMPopu(p,g),bestmem,dim_index);
                tempDecs(:,dim_index) = pop;
            end
            Archive          = [Archive,BestIndividual];            
            BigPopulation{p} = tempDecs;	% record the new position for this small population
            
            if abs(lastBestValA(p)-BestVal) < unChangeThr	% check whether the sub population has converged
                unUpdateNum(p)       = unUpdateNum(p) + 1;
                stopTag(p)           = true;
                ConvergedSolutionSet = [ConvergedSolutionSet,BestIndividual]; 
            end
            lastBestValA(p) = BestVal;
        end
        
        % generate new sub-populations for next stage
        Tag = Global.evaluated > 0.6*Global.evaluation && firstTag;
        if sum(stopTag) == length(stopTag) || Tag 
            firstTag = false;	% The first stage has been over
            
            % evolute the diversity-related variables
            for repPV = 1: 200  % 200 denotes the repeat times for diversity-related variables
                CR      = 0.2;
                F       = 0.5;
                ExiDecs = Archive.decs;
                ExiPV   = ExiDecs(:, PV);
                [N,D]   = size(ExiPV);
                Parent1Dec   = ExiPV;
                Parent2Dec   = ExiPV(randperm(N), :);
                Parent3Dec   = ExiPV(randperm(N), :);
                OffspringDec = Parent1Dec;
                Site = rand(N,D) < CR;
                OffspringDec(Site) = OffspringDec(Site) + F*(Parent2Dec(Site)-Parent3Dec(Site));
                
                Lower = repmat(Global.lower(PV),N,1);	% Lower boundary
                Upper = repmat(Global.upper(PV),N,1);	% Upper boundary
                OffspringDec = max(min(OffspringDec,Upper),Lower);
                
                newDecs        = ExiDecs;
                newDecs(:,PV)  = OffspringDec;
                newPop         = INDIVIDUAL(newDecs);
                
                % environmental selection
                Archive = UpdateArchive(Global.N,[newPop,Archive]);
            end
            % generate new sub-populations
            [BigPopulation,CMAParaMPopu] = GenerateBigPopulation(PV,Groups,Archive);
            
            popN         = length(BigPopulation);
            lastBestValA = 1e+20*ones(1,popN);
            stopTag      = false(1,popN);
            unUpdateNum  = zeros(1,popN);            
            ConvergedSolutionSet = [];
        end
        
        % Obtain the output solutions
        if Global.evaluated >= Global.evaluation
            Decs = zeros(length(BigPopulation),Global.D);
            for bp = 1: length(BigPopulation)
                DecPop      = BigPopulation{bp};
                Decs(bp, :) = DecPop(1,:);
            end
            finalPop = INDIVIDUAL(Decs);
            Archive  = UpdateArchive(Global.N,[finalPop,Archive]);
        end
    end
