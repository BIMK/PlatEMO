function NewParticles = dMOPSO_operator(Global,Particles)
% <operator> <real>
% Particle swarm optimization in dMOPSO
% W --- 0.4 --- The inertia weight in PSO

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    W = Global.ParameterSet(0.4);
    Particles      = Particles([1:end,1:ceil(end/3)*3-end]);
    ParticlesDec   = Particles.decs;
    [N,D]          = size(ParticlesDec);
    ParticlesSpeed = Particles.adds(zeros(N,D));

    %% PSO
    ParticleDec   = ParticlesDec(1:N/3,:);
    ParticleSpeed = ParticlesSpeed(1:N/3,:);
    PBestDec      = ParticlesDec(N/3+1:N/3*2,:);
    GBestDec      = ParticlesDec(N/3*2+1:end,:);
    r1 = repmat(rand(N/3,1),1,D);
    r2 = repmat(rand(N/3,1),1,D);
    NewSpeed = W.*ParticleSpeed + r1.*(PBestDec-ParticleDec) + r2.*(GBestDec-ParticleDec);
    NewDec   = ParticleDec + NewSpeed;
    
    %% Deterministic back
    repair = NewDec < repmat(Global.lower,N/3,1) | NewDec > repmat(Global.upper,N/3,1);
    NewSpeed(repair) = -NewSpeed(repair);

    NewParticles = INDIVIDUAL(NewDec,NewSpeed);
end