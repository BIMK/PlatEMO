function Gbest = WOF_optimiseBySMPSO(GlobalDummy, inputPopulation, evaluations, isDummy)
% This function performs the optimisation by using the SMPSO
% methodology. It is derived from the PlatEMO version of SMPSO.
% Additional functionality has been added and modified to make it
% applicable for normal Individuals as well as the transformed
% weight-variables. 

% ----------------------------------------------------------------------- 
%  WOF_optimiseBySMPSO.m 
%  Copyright (C) 2018 Heiner Zille
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Author of this Code: 
%   Heiner Zille <heiner.zille@ovgu.de>
%
%  This file belongs to the following publications:
%
%  1) Heiner Zille and Sanaz Mostaghim
%     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
%     https://ieeexplore.ieee.org/document/8280974 
% 
%  2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
%     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
%     http://ieeexplore.ieee.org/document/7929324
%  
%  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
%     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
%     http://dl.acm.org/citation.cfm?id=2908979
%
%  Date of publication: 12.10.2018 
%  Last Update: 12.10.2018
% -----------------------------------------------------------------------

% Original copyright disclaimer of the "SMPSO" function of the 
% PlatEMO framework version 1.5: 

    %% Generate random population
    Population       = inputPopulation;
    Pbest            = Population;
    [Gbest,CrowdDis] = WOF_UpdateGbest(Population,GlobalDummy.N);
    
    maximum = currentEvaluations(GlobalDummy, isDummy) + evaluations;

    %% Optimization
    while currentEvaluations(GlobalDummy, isDummy) < maximum
        Population       = WOF_SMPSO_operator(GlobalDummy, [Population,Pbest,Gbest(TournamentSelection(2,GlobalDummy.N,-CrowdDis))], isDummy);
        [Gbest,CrowdDis] = WOF_UpdateGbest([Gbest,Population],GlobalDummy.N);
        Pbest            = WOF_UpdatePbest(Pbest,Population);
    end

end

function e = currentEvaluations(GlobalDummy, isDummy)
    if isDummy == true  
        e = GlobalDummy.Global.evaluated;
    else
        e = GlobalDummy.evaluated;
    end
end

function Pbest = WOF_UpdatePbest(Pbest,Population)
    % Update the local best position of each particle
    replace        = ~all(Population.objs>=Pbest.objs,2);
    Pbest(replace) = Population(replace);
end

function [Gbest,CrowdDis] = WOF_UpdateGbest(Gbest,N)
    % Update the global best set
    Gbest    = Gbest(NDSort(Gbest.objs,1)==1);
    CrowdDis = CrowdingDistance(Gbest.objs);
    [~,rank] = sort(CrowdDis,'descend');
    Gbest    = Gbest(rank(1:min(N,length(Gbest))));
    CrowdDis = CrowdDis(rank(1:min(N,length(Gbest))));
end

function CrowdDis = CrowdingDistance(PopObj)
    % Calculate the crowding distance of each solution in the same front
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end