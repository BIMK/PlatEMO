classdef C10MOP1 < EvoXBenchProblem
% <2024> <multi> <real> <large/none> <expensive/none>
% EvoXBench on CIFAR-10

%------------------------------- Reference --------------------------------
% Z. Lu, R. Cheng, Y. Jin, K. C. Tan, and K. Deb, Neural architecture
% search as multiobjective optimization benchmarks: Problem formulation and
% performance assessment, IEEE Transactions on Evolutionary Computation,
% 2024, 28(2): 323-337.
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            config.name = 'nb101';
            config.args.objs = 'err&params';
            config.args.normalized_objectives = true;
            obj.Setting@EvoXBenchProblem(config);
        end
    end
end

