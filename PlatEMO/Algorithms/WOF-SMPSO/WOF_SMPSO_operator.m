function NewParticles = WOF_SMPSO_operator(Global,Particles, isDummy)
% Particle swarm optimization operator in SMPSO 
% It is derived from the "SMPSO_operator" function of the
% PlatEMO framework and adjusted to fit the WOF algorithm. 

% ----------------------------------------------------------------------- 
%  WOF_SMPSO_operator.m 
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

% Original copyright disclaimer of the "SMPSO_operator" function of the 
% PlatEMO framework version 1.5: 
    
    Particles      = Particles([1:end,1:ceil(end/3)*3-end]);
    ParticlesDec   = Particles.decs;
    [N,D]          = size(ParticlesDec);
    ParticlesSpeed = Particles.adds(zeros(N,D));

    %% PSO
    ParticleDec   = ParticlesDec(1:N/3,:);
    ParticleSpeed = ParticlesSpeed(1:N/3,:);
    PBestDec      = ParticlesDec(N/3+1:N/3*2,:);
    GBestDec      = ParticlesDec(N/3*2+1:end,:);
    W  = repmat(unifrnd(0.1,0.5,N/3,1),1,D);
    r1 = repmat(rand(N/3,1),1,D);
    r2 = repmat(rand(N/3,1),1,D);
    C1 = repmat(unifrnd(1.5,2.5,N/3,1),1,D);
    C2 = repmat(unifrnd(1.5,2.5,N/3,1),1,D);
    NewSpeed = W.*ParticleSpeed + C1.*r1.*(PBestDec-ParticleDec) + C2.*r2.*(GBestDec-ParticleDec);
    phi      = max(4,C1+C2);
    NewSpeed = NewSpeed.*2./abs(2-phi-sqrt(phi.^2-4*phi));
    delta    = repmat((Global.upper-Global.lower)/2,N/3,1);
    NewSpeed = max(min(NewSpeed,delta),-delta);
    NewDec   = ParticleDec + NewSpeed;
    
    %% Deterministic back
    Lower  = repmat(Global.lower,N/3,1);
    Upper  = repmat(Global.upper,N/3,1);
    repair = NewDec < Lower | NewDec > Upper;
    NewSpeed(repair) = 0.001*NewSpeed(repair);
    
    %% Polynomial mutation
    disM  = 20;
    Site1 = repmat(rand(N/3,1)<0.15,1,D);
    Site2 = rand(N/3,D) < 1/D;
    mu    = rand(N/3,D);
    temp  = Site1 & Site2 & mu<=0.5;
    NewDec(temp) = NewDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(NewDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site1 & Site2 & mu>0.5; 
    NewDec(temp) = NewDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-NewDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    if isDummy == true
        NewParticles = [];
        for i = 1:N/3
            NewParticles = [NewParticles, WOF_WeightIndividual(NewDec(i,:),Global,NewSpeed(i,:))];
        end
    else
        NewParticles = INDIVIDUAL(NewDec,NewSpeed);
    end
    
end