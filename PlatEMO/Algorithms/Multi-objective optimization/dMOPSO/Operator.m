function Offspring = Operator(Problem,Particle,Pbest,Gbest,W)
% The particle swarm optimization in dMOPSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin < 5
        W = 0.4;
    end
    ParticleDec = Particle.decs;
    PbestDec    = Pbest.decs;
    GbestDec    = Gbest.decs;
    [N,D]       = size(ParticleDec);
    ParticleVel = Particle.adds(zeros(N,D));

    %% Particle swarm optimization
    r1     = repmat(rand(N,1),1,D);
    r2     = repmat(rand(N,1),1,D);
    OffVel = W.*ParticleVel + r1.*(PbestDec-ParticleDec) + r2.*(GbestDec-ParticleDec);
    OffDec = ParticleDec + OffVel;
    
    %% Deterministic back
    repair = OffDec < repmat(Problem.lower,N,1) | OffDec > repmat(Problem.upper,N,1);
    OffVel(repair) = -OffVel(repair);
    Offspring = Problem.Evaluation(OffDec,OffVel);
end