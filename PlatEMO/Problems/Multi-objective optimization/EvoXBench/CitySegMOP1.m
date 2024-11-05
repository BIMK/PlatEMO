classdef CitySegMOP1 < EvoXBenchProblem
% <2024> <multi> <real> <large/none> <expensive/none>

%------------------------------- Reference --------------------------------
% Z. Lu, R. Cheng, Y. Jin, K. C. Tan, and K. Deb, Neural architecture
% search as multiobjective optimization benchmarks: Problem formulation and
% performance assessment, IEEE Transactions on Evolutionary Computation,
% 2024, 28(2): 323-337.
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            config.name = 'mosegnas';
            config.args.objs = 'err&h1_latency';
            config.args.normalized_objectives = false;
            obj.Setting@EvoXBenchProblem(config);
        end
    end
end