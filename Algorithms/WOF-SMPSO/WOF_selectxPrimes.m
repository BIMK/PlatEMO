function weightIndList = WOF_selectxPrimes(input,amount, method)
% Implements the selection of the x' solutions in WOF-SMPSO. 
% Three methods can be chosen. The first one uses the largest Crowding
% Distance values from the first non-dominated front. The second one
% uses a tournament selection on the population based on
% Pareto-dominance and Crowding Distance. The third option is 
% introduced in publication (1), see above, based on
% reference directions for the first m+1 chosen solutions and selects
% random solutions afterwards.    

% ----------------------------------------------------------------------- 
%  WOF_selectxPrimes.m 
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
    
    inputSize = size(input,2);
    switch method 
        case 1 %largest Crowding Distance from first front
            inFrontNo    = NDSort(input.objs,inf);
            weightIndList = [];
            i = 1;
            if inputSize < amount
                weightIndList = input;
            else
                while size(weightIndList,2) < amount 
                    left = amount - size(weightIndList,2);
                    theFront = inFrontNo == i;
                    newPop = input(theFront);
                    FrontNo    = NDSort(newPop.objs,inf);
                    CrowdDis   = CrowdingDistance(newPop.objs,FrontNo);
                    [~,I] = sort(CrowdDis,'descend');
                    until=min(left,size(newPop,2));
                    weightIndList = [weightIndList,newPop(I(1:until))];
                    i=i+1;
                end
            end
        case 2 %tournament selection by front and CD
            FrontNo    = NDSort(input.objs,inf);
            CrowdDis   = WOF_CrowdingDistance(input.objs,FrontNo);
            weightIndList = input(TournamentSelection(2,amount,FrontNo,-CrowdDis));
        case 3 % first m+1 by reference lines + fill with random
            objValues = input.objs;
            m = size(objValues,2);
            weightIndList = [];
            for i = 1:m
                vec = zeros(1,m);
                vec(1,i) = 1;
                angles = pdist2(vec,objValues,'cosine');
                [minAngle,minIndex] = min(angles);
                weightIndList = [weightIndList,input(minIndex)];
            end
            if size(weightIndList,2) < amount
                vec = ones(1,m);
                angles = pdist2(vec,objValues,'cosine');
                [minAngle,minIndex] = min(angles);
                weightIndList = [weightIndList,input(minIndex)];
            end
            while size(weightIndList,2) < amount
                ind = input(randi([1 inputSize],1,1));
                weightIndList = [weightIndList,ind];
            end
    end
end