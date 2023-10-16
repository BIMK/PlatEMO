function Offspring = Operator(Problem,Particle,PBA,NBA)
% Particle swarm optimization in MO_Ring_PSO_SCD

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    ParticleDec = Particle.decs;
    [N,D]       = size(ParticleDec);
    PbestDec    = zeros(N,D);
    NbestDec    = zeros(N,D);
    for i = 1:N
        PbestDec(i,:) = PBA{i}(1).dec;
        NbestDec(i,:) = NBA{i}(1).dec;        
    end
    ParticleVel = Particle.adds(zeros(N,D));

    %% Particle swarm optimization
    W  = repmat(0.7298,N,D);
    r1 = rand(N,D);
    r2 = rand(N,D);
    C1 = repmat(2.05,N,D);
    C2 = repmat(2.05,N,D);
    OffVel = W.*ParticleVel + C1.*r1.*(PbestDec-ParticleDec) + C2.*r2.*(NbestDec-ParticleDec);
    delta  = repmat((Problem.upper-Problem.lower)/2,N,1);
    OffVel = max(min(OffVel,delta),-delta);
    OffDec = ParticleDec + OffVel;
    Lower  = repmat(Problem.lower,N,1);
    Upper  = repmat(Problem.upper,N,1);
    temp = OffDec<Lower;
    OffDec(temp) = Lower(temp)+0.25*(Upper(temp)-Lower(temp)).*rand(size(OffDec(temp)));
    temp = OffDec>Upper;
    OffDec(temp) = Upper(temp)-0.25*(Upper(temp)-Lower(temp)).*rand(size(OffDec(temp))); 
    Offspring = Problem.Evaluation(OffDec,OffVel);
end
