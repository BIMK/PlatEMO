classdef GLMO < ALGORITHM
% <multi> <real/integer> <large/none>
% Grouped and linked mutation operator algorithm
% optimiser      --- 3 --- The optimisation method. 1 = SMPSO, 2 = NSGA-II, 3 = NSGA-III. Default = NSGA-III
% typeOfGroups   --- 2 --- Grouping method, 1 = linear, 2 = ordered, 3 = random. Default = ordered
% numberOfGroups --- 4 --- The number of varibale Groups. Default = 4  

%------------------------------- Reference --------------------------------
% H. Zille, Large-scale Multi-objective Optimisation: New Approaches and a
% Classification of the State-of-the-Art, PhD Thesis, Otto von Guericke
% University Magdeburg, 2019.
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

    methods
        function main(Algorithm,Problem)
            [optimiser,typeOfGroups,numberOfGroups] = Algorithm.ParameterSet(3,2,4);
            if optimiser == 1
                GLMO_SMPSO(Algorithm,Problem,numberOfGroups,typeOfGroups);
            elseif optimiser == 2
                GLMO_NSGAII(Algorithm,Problem,numberOfGroups,typeOfGroups);
            else
                GLMO_NSGAIII(Algorithm,Problem,numberOfGroups,typeOfGroups);
            end
        end
        function GLMO_SMPSO(Algorithm,Problem,numberOfGroups,typeOfGroups)
            Population       = Problem.Initialization();
            Pbest            = Population;
            [Gbest,CrowdDis] = GLMO_UpdateGbest(Population,Problem.N);
            while Algorithm.NotTerminated(Gbest)
                Population       = GLMO_SMPSOOperator(Problem, Population, Pbest, Gbest(TournamentSelection(2,Problem.N,-CrowdDis)), numberOfGroups, typeOfGroups);
                [Gbest,CrowdDis] = GLMO_UpdateGbest([Gbest,Population],Problem.N);
                Pbest            = GLMO_UpdatePbest(Pbest,Population);
            end
        end
        function GLMO_NSGAII(Algorithm,Problem,numberOfGroups,typeOfGroups)
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = GLMO_NSGAIIEnvironmentalSelection(Population,Problem.N);
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = GLMO_GA(Problem, Population(MatingPool), numberOfGroups, typeOfGroups);
                [Population,FrontNo,CrowdDis] = GLMO_NSGAIIEnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
        function GLMO_NSGAIII(Algorithm,Problem,numberOfGroups,typeOfGroups)
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = GLMO_GA(Problem, Population(MatingPool), numberOfGroups, typeOfGroups);
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = GLMO_NSGAIIIEnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
            end
        end
    end
end