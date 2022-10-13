function [outIndexList,numberOfGroups] = WOF_createGroups(Problem,numberOfGroups,xPrime,method)
% Creates groups of the varibales. Three diffeent methods can be
% chosen. The first one uses linear groups, the second orders variables
% by absolute values, the third is a random grouping. For more
% information about these mechanisms see publications (1) and (3), see below. 
    
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
    
    switch method
        case 1 %linear grouping
            varsPerGroup = floor(Problem.D/numberOfGroups);
            outIndexList = [];
            for i = 1:numberOfGroups-1
               outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
            end
            outIndexList = [outIndexList, ones(1,Problem.D-size(outIndexList,2)).*numberOfGroups];
        case 2 %orderByValueGrouping
            varsPerGroup = floor(Problem.D/numberOfGroups);
            vars = xPrime.dec;
            [~,I] = sort(vars);
            outIndexList = ones(1,Problem.D);
            for i = 1:numberOfGroups-1
               outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
            end
            outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
        case 3 %random Grouping
            varsPerGroup = floor(Problem.D/numberOfGroups);
            outIndexList = [];
            for i = 1:numberOfGroups-1
               outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
            end
            outIndexList = [outIndexList, ones(1,Problem.D-size(outIndexList,2)).*numberOfGroups];
            outIndexList = outIndexList(randperm(length(outIndexList)));
        case 4 %up or down groups
            outIndexList = ones(1,Problem.D);
            xPrimeVars = xPrime.decs;
            xPrimeObjs = xPrime.objs;
            for i = 1 : Problem.D
                newSolVars = xPrime.decs;
                newSolVars(i) = xPrimeVars(i)*1.05;
                newSol = Problem.Evaluation(newSolVars);
                newSolObjs = newSol.objs;
                if newSolObjs(1) < xPrimeObjs(1)
                    outIndexList(i) = 2;
                end
            end
            numberOfGroups = 2;
    end
end