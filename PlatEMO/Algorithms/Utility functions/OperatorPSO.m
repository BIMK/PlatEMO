function Offspring = OperatorPSO(Particle,Pbest,Gbest,W)
%OperatorPSO - The operator of particle swarm optimization.
%
%   Off = OperatorPSO(P,Pbest,Gbest) uses the operator of particle swarm
%   optimization to generate offsprings based on particles P, personal best
%   particles Pbest, and global best particles Gbest. P, Pbest, and Gbest
%   should be arrays of SOLUTION objects, and Off is also an array of
%   SOLUTION objects. Each object of P, Pbest, and Gbest is used to
%   generate one offspring.
%
%   Off = OperatorPSO(P,Pbest,Gbest,W) specifies the parameter of the
%   operator, where W is the inertia weight.
%
%   Example:
%       Off = OperatorPSO(Population,Pbest,Gbest)

%------------------------------- Reference --------------------------------
% C. A. Coello Coello and M. S. Lechuga, MOPSO: A proposal for multiple
% objective particle swarm optimization, Proceedings of the IEEE Congress
% on Evolutionary Computation, 2002, 1051-1056.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin < 4
        W = 0.4;
    end
    ParticleDec = Particle.decs;
    PbestDec    = Pbest.decs;
    GbestDec    = Gbest.decs;
    [N,D]       = size(ParticleDec);
    ParticleVel = Particle.adds(zeros(N,D));

    %% Particle swarm optimization
    r1        = repmat(rand(N,1),1,D);
    r2        = repmat(rand(N,1),1,D);
    OffVel    = W.*ParticleVel + r1.*(PbestDec-ParticleDec) + r2.*(GbestDec-ParticleDec);
    OffDec    = ParticleDec + OffVel;
    Offspring = SOLUTION(OffDec,OffVel);
end