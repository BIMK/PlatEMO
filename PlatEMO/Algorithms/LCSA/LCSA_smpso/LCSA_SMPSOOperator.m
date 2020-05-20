function Offspring = LCSA_SMPSOOperator(Particle,Pbest,Gbest)
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
% This file is derived from its original version containied in the PlatEMO 
% framework. The original copyright disclaimer can be found below. 
% -----------------------------------------------------------------------  

% Particle swarm optimization in SMPSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    ParticleDec = Particle.decs;
    PbestDec    = Pbest.decs;
    GbestDec    = Gbest.decs;
    [N,D]       = size(ParticleDec);
    ParticleVel = Particle.adds(zeros(N,D));
    Global      = GLOBAL.GetObj();

    %% Particle swarm optimization
    W  = repmat(unifrnd(0.1,0.5,N,1),1,D);
    r1 = repmat(rand(N,1),1,D);
    r2 = repmat(rand(N,1),1,D);
    C1 = repmat(unifrnd(1.5,2.5,N,1),1,D);
    C2 = repmat(unifrnd(1.5,2.5,N,1),1,D);
    OffVel = W.*ParticleVel + C1.*r1.*(PbestDec-ParticleDec) + C2.*r2.*(GbestDec-ParticleDec);
    phi    = max(4,C1+C2);
    OffVel = OffVel.*2./abs(2-phi-sqrt(phi.^2-4*phi));
    delta  = repmat((Global.upper-Global.lower)/2,N,1);
    OffVel = max(min(OffVel,delta),-delta);
    OffDec = ParticleDec + OffVel;
    
    %% Deterministic back
    Lower  = repmat(Global.lower,N,1);
    Upper  = repmat(Global.upper,N,1);
    repair = OffDec < Lower | OffDec > Upper;
    OffVel(repair) = 0.001*OffVel(repair);
    OffDec = max(min(OffDec,Upper),Lower);
    
    %% Polynomial mutation
    disM  = 20;
    Site1 = repmat(rand(N,1)<0.15,1,D);
    Site2 = rand(N,D) < 1/D;
    mu    = rand(N,D);
    temp  = Site1 & Site2 & mu<=0.5;
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site1 & Site2 & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = INDIVIDUAL(OffDec,OffVel);
end





