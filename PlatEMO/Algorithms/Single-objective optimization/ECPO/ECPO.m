classdef ECPO < ALGORITHM
% <2021> <single> <real/integer> <large/none> <constrained/none>
% Electric charged particles optimization
% Strategy --- 2 --- Strategy
% nECPI    --- 3 --- Number of ECPs in interaction
% naECP    ---   --- Size of the archive pool

%------------------------------- Reference --------------------------------
% H. R. E. H. Bouchekara. Electric charged particles optimization and its
% application to the optimal design of a circular antenna array.
% Artificial Intelligence Review, 2021, 54: 1767-1802.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by H. R. E. H. Bouchekara (email: bouchekara.houssem@gmail.com)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [Strategy,nECPI,naECP] = Algorithm.ParameterSet(2,3,round(Problem.N/3));

            %% Generate random population
            ECP      = Problem.Initialization();
            [~,rank] = sort(ECP.objs);
            ECP      = ECP(rank);

            if Strategy == 1
                pop_fac = 2*nchoosek(nECPI,2);
            elseif Strategy == 2
                pop_fac = nECPI;
            elseif Strategy == 3
                pop_fac = 2*nchoosek(nECPI,2) + nECPI;
            end

            %% Optimization
            while Algorithm.NotTerminated(ECP)
                Archive   = ECP(1:naECP);
                ECP_NEW   = OperatorECPO(Problem,Strategy,ECP,Archive,pop_fac,nECPI);
                ECP_All   = [Archive, ECP_NEW];
                [~,rank1] = sort(ECP_All.objs);
                ECP       = ECP_All (rank1);
            end
        end
    end
end