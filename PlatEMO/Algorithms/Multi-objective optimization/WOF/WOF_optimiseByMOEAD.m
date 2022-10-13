function Population = WOF_optimiseByMOEAD(GlobalDummy,Population,W,evaluations,isDummy)
% ----------------------------------------------------------------------- 
%  Copyright (C) 2020 Heiner Zille
%
%  This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
%  International License. (CC BY-NC-SA 4.0). To view a copy of this license, 
%  visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or see the 
%  pdf-file "License-CC-BY-NC-SA-4.0.pdf" that came with this code. 
%
%  You are free to: 
%  * Share ? copy and redistribute the material in any medium or format
%  * Adapt ? remix, transform, and build upon the material 
%  Under the following terms:
%  * Attribution ? You must give appropriate credit, provide a link to the 
%     license, and indicate if changes were made. You may do so in any reasonable 
%     manner, but not in any way that suggests the licensor endorses you or your use.
%  * NonCommercial ? You may not use the material for commercial purposes.
%  * ShareAlike ? If you remix, transform, or build upon the material, you must 
%    distribute your contributions under the same license as the original.
%  * No additional restrictions ? You may not apply legal terms or technological 
%    measures that legally restrict others from doing anything the license permits.
% 
%  Author of this Code: 
%   Heiner Zille <heiner.zille@ovgu.de> or <heiner.zille@gmail.com>
%
%  This code is based on the following publications:
% 
%  1) Heiner Zille 
%     "Large-scale Multi-objective Optimisation: New Approaches and a Classification of the State-of-the-Art"  
%     PhD Thesis, Otto von Guericke University Magdeburg, 2019 
%     http://dx.doi.org/10.25673/32063 
% 
%  2) Heiner Zille and Sanaz Mostaghim
%     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
%     https://ieeexplore.ieee.org/document/8280974 
% 
%  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
%     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
%     http://ieeexplore.ieee.org/document/7929324
%  
%  4) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
%     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
%     http://dl.acm.org/citation.cfm?id=2908979
%
%  This file is intended to work with the PlatEMO framework version 2.5. 
%  Date of publication of this code: 06.04.2020 
%  Last Update of this code: 06.04.2020
%  A newer version of this algorithm may be available. Please contact the author 
%  or see http://www.ci.ovgu.de/Research/Codes.html. 
%
% The files may have been modified in Feb 2021 by the authors of the Platemo framework to work with the Platemo 3.0 release. 
% ----------------------------------------------------------------------- 
% This file is derived from its original version containied in the PlatEMO 
% framework.  
% -----------------------------------------------------------------------

% MOEA/D

    maximum = currentEvaluations(GlobalDummy, isDummy) + evaluations;
    
    T = max(ceil(GlobalDummy.N/10),2);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    
    %% Associate each subproblem with one solution
    % The ideal point
    Z = min(Population.objs,[],1);
    % The value of each solution on each subproblem (modified Tchebycheff approach)
    g = zeros(GlobalDummy.N);
    for i = 1 : GlobalDummy.N
        g(i,:) = max(repmat(abs(Population(i).obj-Z),GlobalDummy.N,1)./W,[],2)';
    end
    [~,rank] = sort(g,2);
    % The index of solution which each subproblem associated with
    associate = zeros(1,GlobalDummy.N);
    for i = 1 : GlobalDummy.N
        x = find(~associate(rank(i,:)),1);
        associate(rank(i,x)) = i;
    end
    Population = Population(associate);
    
    %% Optimization
    while currentEvaluations(GlobalDummy, isDummy) < maximum
        % For each solution
        for i = 1 : GlobalDummy.N
            % Choose the parents
            P = B(i,randperm(size(B,2)));

            % Generate an offspring
            if isDummy == true
                NewDec = WOF_GAhalf(GlobalDummy, Population(P(1:2)));
                Offspring = WOF_WeightIndividual(NewDec,GlobalDummy);
        
            else
                Offspring = OperatorGAhalf(GlobalDummy,Population(P(1:2)));
            end
            
            % Update the ideal point
            Z = min(Z,Offspring.obj);

            type = 1;
            % Update the neighbours
            switch type
                case 1
                    % PBI approach
                    normW   = sqrt(sum(W(P,:).^2,2));
                    normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                    normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                    CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                    CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                    g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                case 2
                    % Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                case 3
                    % Tchebycheff approach with normalization
                    Zmax  = max(Population.objs,[],1);
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                case 4
                    % Modified Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
            end
            Population(P(g_old>=g_new)) = Offspring;
        end
    end
end

function e = currentEvaluations(GlobalDummy, isDummy)
    if isDummy == true  
        e = GlobalDummy.Global.FE;
    else
        e = GlobalDummy.FE;
    end
end