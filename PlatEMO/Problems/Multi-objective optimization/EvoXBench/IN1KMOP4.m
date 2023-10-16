classdef IN1KMOP4 < EvoXBenchProblem
% <multi> <real> <large/none> <expensive/none>
% EvoXBench on ImageNet 1K

%------------------------------- Reference --------------------------------
% Z. Lu, R. Cheng, Y. Jin, K. C. Tan, and K. Deb, Neural architecture
% search as multiobjective optimization benchmarks: Problem formulation and
% performance assessment, IEEE Transactions on Evolutionary Computation,
% 2023.
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            config.name = 'transformer';
            config.args.objs = 'err&params';
            config.args.normalized_objectives = false;
            obj.Setting@EvoXBenchProblem(config);
        end
    end
end
