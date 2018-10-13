function NewParticles = MMOPSO_operator(Global,Particles)
% <operator> <real>
% Particle swarm optimization in MMOPSO

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Particles      = Particles([1:end,1:ceil(end/3)*3-end]);
    ParticlesDec   = Particles.decs;
    [N,D]          = size(ParticlesDec);
    ParticlesSpeed = Particles.adds(zeros(N,D));

    %% PSO
    ParticleDec   = ParticlesDec(1:N/3,:);
    ParticleSpeed = ParticlesSpeed(1:N/3,:);
    PBestDec      = ParticlesDec(N/3+1:N/3*2,:);
    GBestDec      = ParticlesDec(N/3*2+1:end,:);
    W  = repmat(rand(N/3,1)*0.4+0.1,1,D);
    r1 = repmat(rand(N/3,1),1,D);
    r2 = repmat(rand(N/3,1),1,D);
    C1 = repmat(rand(N/3,1)*0.5+1.5,1,D);
    C2 = repmat(rand(N/3,1)*0.5+1.5,1,D);
    NewSpeed        = W.*ParticleSpeed;
    temp            = repmat(rand(N/3,1)<0.7,1,D);
    NewSpeed(temp)  = NewSpeed(temp) + C1(temp).*r1(temp).*(PBestDec(temp)-ParticleDec(temp));
    NewSpeed(~temp) = NewSpeed(~temp) + C2(~temp).*r2(~temp).*(GBestDec(~temp)-ParticleDec(~temp));
    NewDec          = ParticleDec + NewSpeed;
    
    NewParticles = INDIVIDUAL(NewDec,NewSpeed);
end