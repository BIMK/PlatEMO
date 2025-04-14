function Offspring = OperatorAESSPSO(Problem,Particle,Pbest,Gbest,Beta, Gamma)
%OperatorPSO - The operator of particle swarm optimization.
%
%   Off = OperatorAESSPSO(Pro, P, Pbest, Gbest, Beta, Gamma) uses the operator of particle
%   swarm optimization to generate offsprings for problem Pro based on
%   particles P, personal best particles Pbest, and global best particles
%   Gbest. P, Pbest, and Gbest should be arrays of SOLUTION objects, and
%   Off is also an array of SOLUTION objects. Each object of P, Pbest, and
%   Gbest is used to generate one offspring. 
%   
%   Beta and Gamma are control parameters that influence the exploration and
%   exploitation behavior of the algorithm. Beta adjusts the weight of the
%   personal best component, while Gamma scales the influence of the global
%   best component. These parameters allow for adaptive tuning of the
%   algorithm's search dynamics.
%
%   Off = OperatorAESSPSO(Pro,P,Pbest,Gbest,Beta, Gamma) specifies the parameter of the
%   operator.
%
%   Example:
%       Off = OperatorAESSPSO(Problem,Population,Pbest,Gbest,Beta, Gamma)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Mehdi Alimohammadi (alimohammadimd@chmail.ir)

    %% Parameter setting
    if nargin < 5
        Beta  = 2.05;
        Gamma = 2.05;
    end
    ParticleDec = Particle.decs;
    PbestDec    = Pbest.decs;
    GbestDec    = Gbest.decs;
    [N,D]       = size(ParticleDec);
    ParticleVel = Particle.adds(zeros(N,D));

    %% Particle swarm optimization
    r1 = Beta*rand(N,1);
    r2 = Gamma*rand(N,1);
    W  = w_adapter(N,r1,r2);
    
    OffVel    = W.*ParticleVel + r1.*(PbestDec-ParticleDec) + r2.*(GbestDec-ParticleDec);
    OffDec    = ParticleDec + OffVel;
    Offspring = Problem.Evaluation(OffDec,OffVel);
end

function w = w_adapter(n,Beta,Gamma)
    W_min = -10;
    W_max = 10;
    diff  = 1;
    Prec  = 0.0001;
    while diff > Prec
        Omega = (W_min+W_max)/2;
        Omega_left = Omega-0.0000001;
        E   = ctrleigall(n,Beta,Gamma,Omega);
        Min = min(E);
        E   = ctrleigall(n,Beta,Gamma,Omega_left);
        Min_left = min(E);
        if Min > Min_left
            diff  = Omega-W_min;
            W_min = Omega;
        else
            diff  = W_max-Omega;
            W_max = Omega;
        end
    end
    w = (W_min+W_max)/2;
end

function E = ctrleigall(n,Beta,Gamma,Omega)
    B = [diag(Beta) Gamma
         diag(Beta) Gamma];
    A = [diag(1-Beta-Gamma) Omega*eye(n)
         diag(-Beta-Gamma) Omega*eye(n)];
    Q  = B;
    Q1 = A*B;
    Q  = [Q Q1];
    E  = svd(Q);
end