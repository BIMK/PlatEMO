function Offspring = Operator(Particle,Pbest,Gbest)
% Particle swarm optimization in MMOPSO

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

    %% Particle swarm optimization
    W  = repmat(unifrnd(0.1,0.5,N,1),1,D);
    r1 = repmat(rand(N,1),1,D);
    r2 = repmat(rand(N,1),1,D);
    C1 = repmat(unifrnd(1.5,2,N,1),1,D);
    C2 = repmat(unifrnd(1.5,2,N,1),1,D);
    OffVel        = W.*ParticleVel;
    temp          = repmat(rand(N,1)<0.7,1,D);
    OffVel(temp)  = OffVel(temp) + C1(temp).*r1(temp).*(PbestDec(temp)-ParticleDec(temp));
    OffVel(~temp) = OffVel(~temp) + C2(~temp).*r2(~temp).*(GbestDec(~temp)-ParticleDec(~temp));
    OffDec        = ParticleDec + OffVel;
    Offspring     = INDIVIDUAL(OffDec,OffVel);
end