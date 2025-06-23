classdef CoDEA < ALGORITHM
% <2022> <multi/many> <real/integer/label/binary/permutation>
% collaborative decomposition-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% Yu Wu, Jianle Wei, Weiqin Ying, Yanqi Lan, Zhen Cui, Zhenyu Wang,
% A collaborative decomposition-based evolutionary algorithm integrating 
% normal and penalty-based boundary intersection methods for many-objective 
% optimization, Information Sciences,Volume 616,2022,Pages 505-525
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate the reference points and random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population   = Problem.Initialization();
            [z,znad]     = deal(min(Population.objs),max(Population.objs));
            H1 = 1;
            while nchoosek(H1+Problem.M,Problem.M-1) <= Problem.N
                H1 = H1 + 1;
            end
            ILid = nchoosek(H1+Problem.M-1,Problem.M-1);
            r = [];
            for i = 1 : ILid
                deta = min(W(i,:));
                beta = (1 - max(W(i,:)));
                r(i) = ((2*(1-deta*Problem.M)+(beta/(0.5)))/2);%Ðý×ªÏµÊý
            end
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = randi(Problem.N,1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,z,znad] = EnvironmentalSelection([Population,Offspring],W,Problem.N,z,znad,ILid,r );
            end
        end
    end
end