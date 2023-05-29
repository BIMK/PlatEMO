classdef C10MOP5 < EvoXBenchProblem
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
            config.name = 'nb201';
            config.args.fidelity = 200;
            config.args.objs = ['err&params&flops&edgegpu_latency&' ...
                'edgegpu_energy'];
            config.args.dataset = 'cifar10';
            config.args.normalized_objectives = true;
            obj.Setting@EvoXBenchProblem(config);
        end
    end
end