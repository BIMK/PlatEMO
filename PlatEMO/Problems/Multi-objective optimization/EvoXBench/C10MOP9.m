classdef C10MOP9 < EvoXBenchProblem
% <multi> <real> <large/none> <expensive/none>
% EvoXBench on CIFAR-10

%------------------------------- Reference --------------------------------
% Z. Lu, R. Cheng, Y. Jin, K. C. Tan, and K. Deb, Neural architecture
% search as multiobjective optimization benchmarks: Problem formulation and
% performance assessment, IEEE Transactions on Evolutionary Computation,
% 2023.
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            config.name = 'darts';
            config.args.objs = 'err&params&flops';
            config.args.normalized_objectives = false;
            obj.Setting@EvoXBenchProblem(config);
        end
    end
end