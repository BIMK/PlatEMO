function outIndexList = WOF_createGroups(numberOfGroups,xPrime,numberOfVariables, method)
% Creates groups of the varibales. Three diffeent methods can be
% chosen. The first one uses linear groups, the second orders variables
% by absolute values, the third is a random grouping. For more
% information about these mechanisms see publication (2), see above. 
    
% ----------------------------------------------------------------------- 
%  WOF_createGroups.m 
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
    
    switch method
        case 1 %linear grouping
            varsPerGroup = floor(numberOfVariables/numberOfGroups);
            outIndexList = [];
            for i = 1:numberOfGroups-1
               outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
            end
            outIndexList = [outIndexList, ones(1,numberOfVariables-size(outIndexList,2)).*numberOfGroups];
        case 2 %orderByValueGrouping
            varsPerGroup = floor(numberOfVariables/numberOfGroups);
            vars = xPrime.dec;
            [~,I] = sort(vars);
            outIndexList = ones(1,numberOfVariables);
            for i = 1:numberOfGroups-1
               outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
            end
            outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
        case 3 %random Grouping
            varsPerGroup = floor(numberOfVariables/numberOfGroups);
            outIndexList = [];
            for i = 1:numberOfGroups-1
               outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
            end
            outIndexList = [outIndexList, ones(1,numberOfVariables-size(outIndexList,2)).*numberOfGroups];
            outIndexList = outIndexList(randperm(length(outIndexList)));
    end
end