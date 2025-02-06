classdef DKCA < ALGORITHM
% <2024> <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Dynamic knowledge-guided coevolutionary algorithm
% t --- 0.3 --- Proportion of the decision variable to the original variable after dimensionality reduction
% k ---   4 --- Successive k generations of variables change their corresponding scores

%------------------------------- Reference --------------------------------
% Y. Li, X. Feng, and H. Yu. A dynamic knowledge-guided coevolutionary
% algorithm for large-scale sparse multiobjective optimization problems.
% IEEE Transactions on Systems, Man, and Cybernetics: Systems, 2024,
% 54(11): 7054-7064.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            N2    = 4;
            [t,k] = Algorithm.ParameterSet(0.3,4);
            T     = t*Problem.D;
            non_zero   = zeros(1,Problem.N); 
            generation = 1;
            TDec       = [];
            TMask      = [];
            TempPop    = [];
            TMask_d1   = [];
            TempPop_d1 = [];
            Fitness    = zeros(1,Problem.D); 
            REAL       = ~strcmp(Problem.encoding,'binary');  

            %% Reduction and calculate the scores of variables
            for i = 1 : N2
                if REAL
                    Dec  = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                    Dec2 = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                else
                    Dec  = ones(Problem.D,Problem.D);
                    Dec2 = ones(Problem.D,Problem.D);
                end
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            dim         = DimSelection(Problem,Mask,Dec,Dec2,Problem.D,T);  
            dim(dim==0) = [];
            d1          = size(dim,2);  
            init_time   = 4;

            %% Generate initial RP
            for i = 1 : init_time
                Mask_dim = eye(d1);   %d1*d1
                Mask_d1  = zeros(d1,Problem.D);   %d1*D
                for j = 1 : d1
                     Mask_d1(j,dim(find(Mask_dim(j,:)))) = 1;
                end
                if REAL
                    Dec_d1 = unifrnd(repmat(Problem.lower,d1,1),repmat(Problem.upper,d1,1));
                else
                    Dec_d1 = ones(d1,Problem.D);
                end
                Population_d1 = Problem.Evaluation(Dec_d1.*Mask_d1);  
                TMask_d1      = [TMask_d1;Mask_dim];
                TempPop_d1    = [TempPop_d1,Population_d1]; 
            end
            if init_time*d1 < Problem.N
                num_e = init_time*d1;
            else
                num_e = Problem.N;
            end
            [Population_d1,Mask_dim,FrontNo_d1,CrowdDis_d1] = EnvironmentalSelection_withoutDec([Population_d1,TempPop_d1],[Mask_dim;TMask_d1],num_e);

            %% Generate initial OP
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1)); 
            else
                Dec = ones(Problem.N,Problem.D);
            end
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);
           
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool       = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness,REAL);

                MatingPool_d1 = TournamentSelection(2,2*Problem.N,FrontNo_d1,-CrowdDis_d1);
                [~,Mask_dim_off] = Operator(Problem,Dec(MatingPool_d1,:),Mask_dim(MatingPool_d1,:),Fitness,REAL);  

                for j = 1 : Problem.N
                    OffMask(j,dim(find(Mask_dim_off(j,:)))) = 1;
                end
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population_d1,Mask_dim,FrontNo_d1,CrowdDis_d1] = EnvironmentalSelection_withoutDec([Population_d1,Offspring],[Mask_dim;Mask_dim_off],Problem.N);
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                for i = 1 : size(Mask_dim,1)
                    non_zero_temp = find(Mask_dim(i,:));  
                    non_zero(1,i) = size(non_zero_temp,2);  
                end
                non_zero_num = mode(non_zero(1,:),2); 
                if generation > 1   
                    if num_temp == non_zero_num   
                        same_sparsity = same_sparsity+1;
                    else
                        same_sparsity = 1;
                    end
                else
                    same_sparsity = 1;
                end
                num_temp   = non_zero_num;  
                generation = generation+1;
                if same_sparsity > k  
                    dim_base = DimBase(Mask_dim,dim,num_temp);
                    for loc = dim_base  
                        Fitness(loc) = Fitness(loc) - ceil(1/N2*Fitness(loc));
                    end
                end
            end
        end
    end
end