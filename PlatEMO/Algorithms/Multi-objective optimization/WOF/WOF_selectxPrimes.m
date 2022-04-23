function weightIndList = WOF_selectxPrimes(input,amount, method)
% Implements the selection of the x' solutions in WOF-SMPSO. 
% Three methods can be chosen. The first one uses the largest Crowding
% Distance values from the first non-dominated front. The second one
% uses a tournament selection on the population based on
% Pareto-dominance and Crowding Distance. The third option is 
% introduced in publication (2), see above, based on
% reference directions for the first m+1 chosen solutions and selects
% random solutions afterwards.    

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
                    CrowdDis   = WOF_CrowdingDistance(newPop.objs,FrontNo);
                    [~,I] = sort(CrowdDis,'descend');
                    until=min(left,size(newPop,2));
                    weightIndList = [weightIndList,newPop(I(1:until))];
                    i=i+1;
                end
            end
        case 2 %tournament selection by front and CD
            FrontNo    = NDSort(input.objs,inf);
            CrowdDis   = CrowdingDistance(input.objs,FrontNo);
            weightIndList = input(TournamentSelection(2,amount,FrontNo,-CrowdDis));
        case 3 % first m+1 by reference lines + fill with random
            objValues = input.objs;
            m = size(objValues,2);
            weightIndList = [];
            for i = 1:m
                vec = zeros(1,m);
                vec(1,i) = 1;
                angles = pdist2(vec,real(objValues),'cosine');
                [minAngle,minIndex] = min(angles);
                weightIndList = [weightIndList,input(minIndex)];
            end
            if size(weightIndList,2) < amount
                vec = ones(1,m);
                angles = pdist2(vec,real(objValues),'cosine');
                [minAngle,minIndex] = min(angles);
                weightIndList = [weightIndList,input(minIndex)];
            end
            while size(weightIndList,2) < amount
                ind = input(randi([1 inputSize],1,1));
                weightIndList = [weightIndList,ind];
            end
    end
end