classdef CitySegMOP10 < EvoXBenchProblem
% <multi> <real> <large/none> <expensive/none>

%------------------------------- Reference --------------------------------
% Z. Lu, R. Cheng, Y. Jin, K. C. Tan, and K. Deb, Neural architecture
% search as multiobjective optimization benchmarks: Problem formulation and
% performance assessment, IEEE Transactions on Evolutionary Computation,
% 2023.
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            config.name = 'mosegnas';
            config.args.objs = 'err&h2_latency&h2_energy_consumption&flops&params';
            config.args.normalized_objectives = false;
            obj.Setting@EvoXBenchProblem(config);
        end
    end
end

