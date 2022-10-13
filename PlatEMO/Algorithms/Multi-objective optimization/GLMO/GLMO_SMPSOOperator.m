function Offspring = GLMO_SMPSOOperator(Problem, Particle, Pbest, Gbest, numberOfGroups, typeOfGroups)
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
% This file is derived from its original version containied in the PlatEMO 
% framework. The original copyright disclaimer can be found below. 
% ----------------------------------------------------------------------- 

% Particle swarm optimization in SMPSO

    %% Parameter setting
    ParticleDec = Particle.decs;
    PbestDec    = Pbest.decs;
    GbestDec    = Gbest.decs;
    [N,D]       = size(ParticleDec);
    ParticleVel = Particle.adds(zeros(N,D));

    %% Particle swarm optimization
    W  = repmat(unifrnd(0.1,0.5,N,1),1,D);
    r1 = repmat(rand(N,1),1,D);
    r2 = repmat(rand(N,1),1,D);
    C1 = repmat(unifrnd(1.5,2.5,N,1),1,D);
    C2 = repmat(unifrnd(1.5,2.5,N,1),1,D);
    OffVel = W.*ParticleVel + C1.*r1.*(PbestDec-ParticleDec) + C2.*r2.*(GbestDec-ParticleDec);
    phi    = max(4,C1+C2);
    OffVel = OffVel.*2./abs(2-phi-sqrt(phi.^2-4*phi));
    delta  = repmat((Problem.upper-Problem.lower)/2,N,1);
    OffVel = max(min(OffVel,delta),-delta);
    OffDec = ParticleDec + OffVel;
    
    %% Deterministic back
    Lower  = repmat(Problem.lower,N,1);
    Upper  = repmat(Problem.upper,N,1);
    repair = OffDec < Lower | OffDec > Upper;
    OffVel(repair) = 0.001*OffVel(repair);
    OffDec = max(min(OffDec,Upper),Lower);
    
    %% Polynomial mutation
    disM  = 20;
    Site1 = repmat(rand(N,1)<0.15,1,D);

    [outIndexList,~] = GLMO_createGroups(numberOfGroups,OffDec,D,typeOfGroups); % 3 = random groups
    chosengroups = randi(numberOfGroups,size(outIndexList,1),1);
    Site2 = outIndexList == chosengroups;
    mu    = rand(N,1);
    mu    = repmat(mu,1,D);       
    
    temp  = Site1 & Site2 & mu<=0.5;
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site1 & Site2 & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(OffDec,OffVel);
end