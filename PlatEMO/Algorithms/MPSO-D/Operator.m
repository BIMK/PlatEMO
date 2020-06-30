function Offspring = Operator(Particle,Pbest,Gbest,Parameter)
% <operator> <real>
% Particle swarm optimization in MPSO/D
% c1   ---   2 --- Parameter in updating particle's velocity
% c2   ---   2 --- Parameter in updating particle's velocity
% CR   --- 0.5 --- Parameter CR in differental evolution
% F    --- 0.5 --- Parameter F in differental evolution
% proM ---   1 --- The expectation of number of bits doing mutation 
% disM ---  20 --- The distribution index of polynomial mutation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 3
        [c1,c2,CR,F,proM,disM] = deal(Parameter{:});
    else
        [c1,c2,CR,F,proM,disM] = deal(2,2,0.5,0.5,1,20);
    end
    ParticleDec = Particle.decs;
    PbestDec    = Pbest.decs;
    GbestDec    = Gbest.decs;
    [N,D]       = size(ParticleDec);
    ParticleVel = Particle.adds(zeros(N,D));
    Global      = GLOBAL.GetObj();   
    
    %% Particle swarm optimization
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    DoPSO = repmat(rand(N,1)<0.5,1,D);
    W     = 0.9 - Global.evaluated./Global.evaluation*0.8;
    r1    = repmat(rand(N,1),1,D);
    r2    = repmat(rand(N,1),1,D);
    OffVel        = ParticleVel;
    OffDec        = ParticleDec;
    OffVel(DoPSO) = W.*ParticleVel(DoPSO) + c1.*r1(DoPSO).*(PbestDec(DoPSO)-ParticleDec(DoPSO)) + c2.*r2(DoPSO).*(GbestDec(DoPSO)-ParticleDec(DoPSO));
    OffDec(DoPSO) = ParticleDec(DoPSO) + OffVel(DoPSO);
    % Set the infeasible decision variables to the value of their parents
    Invalid         = OffDec < Lower | OffDec > Upper;
    OffDec(Invalid) = ParticleDec(Invalid);
    
    %% DE
    Site = ~DoPSO & rand(N,D)<CR;
    OffDec(Site) = ParticleDec(Site) + F.*(GbestDec(Site)-PbestDec(Site));
    % Set the infeasible decision variables to boundary values
    OffDec = max(min(OffDec,Upper),Lower);

    %% Polynomial mutation
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = INDIVIDUAL(OffDec,OffVel);
end