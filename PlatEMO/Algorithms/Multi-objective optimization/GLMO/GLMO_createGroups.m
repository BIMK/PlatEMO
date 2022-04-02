function [outIndexArray,numberOfGroupsArray] = GLMO_createGroups(numberOfGroups, xPrime, numberOfVariables, method)
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
%  2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Mutation Operators Based on Variable Grouping for Multi-objective Large-scale Optimization"
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Athens, Greece, December 2016
%     https://ieeexplore.ieee.org/document/7850214 
%
%  This file is intended to work with the PlatEMO framework version 2.5. 
%  Date of publication of this code: 06.04.2020 
%  Last Update of this code: 06.04.2020
%  A newer version of this algorithm may be available. Please contact the author 
%  or see http://www.ci.ovgu.de/Research/Codes.html. 
%
% The files may have been modified in Feb 2021 by the authors of the Platemo framework to work with the Platemo 3.0 release. 
% -----------------------------------------------------------------------

% works for variables only, no INDIVIDUALS needed. 
% works also for arrays of solution variables. 

outIndexArray = [];
numberOfGroupsArray = [];

noOfSolutions = size(xPrime,1);
for sol = 1:noOfSolutions
    
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
            vars = xPrime(sol,:);
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
    outIndexArray = [outIndexArray;outIndexList];
    numberOfGroupsArray = [numberOfGroupsArray;numberOfGroups];
    
end
end