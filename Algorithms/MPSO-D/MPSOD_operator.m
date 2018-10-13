function NewParticles = MPSOD_operator(Global,Particles)
% <operator> <real>
% Particle swarm optimization in MPSO/D
% c1   ---   2 --- Parameter in updating particle's velocity
% c2   ---   2 --- Parameter in updating particle's velocity
% CR   --- 0.5 --- Parameter CR in differental evolution
% F    --- 0.5 --- Parameter F in differental evolution
% proM ---   1 --- The expectation of number of bits doing mutation 
% disM ---  20 --- The distribution index of polynomial mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [c1,c2,CR,F,proM,disM] = Global.ParameterSet(2,2,0.5,0.5,1,20);
    Particles      = Particles([1:end,1:ceil(end/3)*3-end]);
    ParticlesDec   = Particles.decs;
    [N,D]          = size(ParticlesDec);
    ParticlesSpeed = Particles.adds(zeros(N,D));

    %% Initialize the parents and offsprings
    % Parents and their pbest and gbest
    ParticleDec   = ParticlesDec(1:N/3,:);
    ParticleSpeed = ParticlesSpeed(1:N/3,:);
    PBestDec      = ParticlesDec(N/3+1:N/3*2,:);
    GBestDec      = ParticlesDec(N/3*2+1:end,:);
    % Initialize the offsprings
    NewSpeed = ParticleSpeed;
    NewDec   = ParticleDec;
    
    %% PSO
    Lower = repmat(Global.lower,N/3,1);
    Upper = repmat(Global.upper,N/3,1);
    DoPSO = repmat(rand(N/3,1)<0.5,1,D);
    W     = 0.9 - Global.evaluated./Global.evaluation*0.8;
    r1    = repmat(rand(N/3,1),1,D);
    r2    = repmat(rand(N/3,1),1,D);
    NewSpeed(DoPSO) = W.*ParticleSpeed(DoPSO) + c1.*r1(DoPSO).*(PBestDec(DoPSO)-ParticleDec(DoPSO)) + c2.*r2(DoPSO).*(GBestDec(DoPSO)-ParticleDec(DoPSO));
    NewDec(DoPSO)   = ParticleDec(DoPSO) + NewSpeed(DoPSO);
    % Set the infeasible decision variables to the value of their parents
    Invalid         = NewDec < Lower | NewDec > Upper;
    NewDec(Invalid) = ParticleDec(Invalid);
    
    %% DE
    Site = ~DoPSO & rand(N/3,D)<CR;
    NewDec(Site) = ParticleDec(Site) + F.*(GBestDec(Site)-PBestDec(Site));
    % Set the infeasible decision variables to boundary values
    NewDec = max(min(NewDec,Upper),Lower);

    %% Polynomial mutation
    Site  = rand(N/3,D) < proM/D;
    mu    = rand(N/3,D);
    temp  = Site & mu<=0.5;
    NewDec(temp) = NewDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(NewDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    NewDec(temp) = NewDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-NewDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    NewParticles = INDIVIDUAL(NewDec,NewSpeed);
end