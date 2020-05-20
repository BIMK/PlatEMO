function LCSA_NSGAIII(Global)
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
%     "Linear Search Mechanism for Multi- and Many-Objective Optimisation"
%     10th International Conference on Evolutionary Multi-Criterion Optimization (EMO 2019), 
%        Lecture Notes in Computer Science, vol 11411. 
%        Deb K. et al. (eds), Springer, Cham, East Lansing, Michigan, USA, March 2019  
%     https://doi.org/10.1007/978-3-030-12598-1_32.
%
%  This file is intended to work with the PlatEMO framework version 2.5. 
%  Date of publication of this code: 06.04.2020 
%  Last Update of this code: 06.04.2020
%  A newer version of this algorithm may be available. Please contact the author 
%  or see http://www.ci.ovgu.de/Research/Codes.html. 
% ----------------------------------------------------------------------- 

    %% Generate the reference points and random population
    [Z,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    Zmin         = min(Population(all(Population.cons<=0,2)).objs,[],1);
    
    iter = 1; 
    xintervall = LCSA_getNumberOfGenerationsForNormalOptimisation(Global.evaluation, Global.N);

    %% Optimization
    while Global.NotTermination(Population)
        
        MatingPool = TournamentSelection(2,Global.N,sum(max(0,Population.cons),2));
        Offspring  = GA(Population(MatingPool));
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        Population = LCSA_NSGAIIIEnvironmentalSelection([Population,Offspring],Global.N,Z,Zmin);
        
        iter = iter + 1;
        if mod(iter,xintervall) == 0
            Global.NotTermination(Population);
            noOfXOptGenerations = LCSA_getNumberOfGenerationsForCoefficientOptimisation(Global.evaluation, Global.N);
            decs                = Population.decs;
            xlower              = ones(1,Global.N)*-10; %-10
            xupper              = ones(1,Global.N)*10; %10
            newPopSize          = Global.N;
            xDecs               = repmat(xlower,newPopSize,1) + rand(newPopSize,Global.N).*repmat(xupper-xlower,newPopSize,1); %newPopSize solutions with N vars each.

            xPop                = INDIVIDUAL(transpose(transpose(decs)*transpose(xDecs)), xDecs);
            
            Zmin                = min([Zmin;xPop(all(xPop.cons<=0,2)).objs],[],1);
            Population          = LCSA_NSGAIIIEnvironmentalSelection([Population,xPop],Global.N,Z,Zmin);
            
            
            
            
            for temp = 1:noOfXOptGenerations 
                
                Global.NotTermination(Population);
                xMatingPool = TournamentSelection(2,Global.N,sum(max(0,xPop.cons),2));
                xDecs       = LCSA_GA(xPop(xMatingPool).adds,{1,20,1,20},xlower,xupper);
                newSols     = INDIVIDUAL(transpose(transpose(decs)*transpose(xDecs)), xDecs);

                Zmin        = min([Zmin;newSols(all(newSols.cons<=0,2)).objs],[],1);
                xPop        = LCSA_NSGAIIIEnvironmentalSelection([xPop,newSols],Global.N,Z,Zmin);
                Population  = LCSA_NSGAIIIEnvironmentalSelection([Population,xPop],Global.N,Z,Zmin);
                
            end 
             
            
        end

        
    end
end