classdef MultiObjectiveEGO < ALGORITHM
% <multi> <real/integer> <constrained/none> <expensive>
% Multi-objective efficient global optimization
% alpha --- 0.7 --- portion of samples for Kriging construction
% num_k ---   5 --- number of infill points per iteration
% H     ---  21 --- number of reference directions

%------------------------------- Reference --------------------------------
% R. Hussein, K. Deb, A Generative Kriging Surrogate Model for Constrained 
% and Unconstrained Multi-objective Optimization, in: Proc. Genet. Evol. 
% Comput. Conf. 2016, Denver, 2016, 573-580. 
%--------------------------------------------------------------------------

% This function is written by Youwei He (email: 1554748356@qq.com)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [alpha,num_k,H] = Algorithm.ParameterSet(0.7,5,21);
            % parameter for AASF in equation (10)
            rho = 1e-3;
            %% Generate the initial design points
            % number of design variables
            D = Problem.D;
            % number of initial design points 11*D-1
            N = Problem.N;
            %% step1-1: generate initial design points using Latin Hypercube sampling
            PopDec = repmat(Problem.upper-Problem.lower,N,1).*UniformPoint(N,D,'Latin') ...
                +repmat(Problem.lower,N,1);
            %% step1-2: evaluate initial design points
            Population = Problem.Evaluation(PopDec);

            %% step 2: Generate the reference direction set
            [R,N_R] = UniformPoint(H,Problem.M);
            R = R./sqrt(sum(R.^2,2));
            %% Optimization
            while Algorithm.NotTerminated(Population)
                %% step 2 diversity_preserver procedure: neighborhood approach
                i_direction = linspace(1,N_R,N_R);
                for i = i_direction
                    for j = 1 : num_k
                        Algorithm.NotTerminated(Population);
                        PopDecT = Population.decs;
                        PopObjT = Population.objs;
                        PopConT = Population.cons;% cons>=0: meet constraints
                        num_con = size(PopConT,2);
                        N_Pop=length(PopDecT(:,1));
                        % determine whether the problem is constrained or not
                        if all(all(PopConT==0)) && N_Pop==1
                            constrained=1;
                        else
                            constrained=0;
                        end
                        %% step 3: Points_Selector procedure
                        normW = sqrt(sum(R(i,:).^2,2));
                        normP = sqrt(sum(PopObjT.^2,2));
                        CosineP = sum((PopObjT).*R(i,:),2)./normP./normW;
                        % orthogonal distance of each point to the given reference direction
                        distB = normP.*sqrt(1-CosineP.^2);
                        [~,indx] = sort(distB);
                        N = ceil(alpha*N_Pop);
                        PopObj = PopObjT(indx(1:N),:);
                        PopDec = PopDecT(indx(1:N),:);
                        PopCon = PopConT(indx(1:N),:);
                        %% step 4-1
                        index = sum(PopCon >= 0, 2) == num_con;% feasible solutions
                        PopCon(index,:)=0;
                        PopCon(~index,:)=-PopCon(~index,:);     
                        PopObjScaled = (PopObj-repmat(min(PopObj),N,1))./(repmat(max(PopObj),N,1)-repmat(min(PopObj),N,1));  
                        if constrained
                            PopConScaled = (PopCon-repmat(min(PopCon),N,1))./(repmat(max(PopCon),N,1)-repmat(min(PopCon),N,1));
                        end
                        %% step 4-2
                        PopSmetric = zeros(N,1);
                        PopSmetricT = PopObjScaled(index,:)./repmat(R(i,:),sum(index),1);
                        PopSmetric(index,:)=max(PopSmetricT,[],2)+rho*sum(PopSmetricT,2);
                        if constrained
                            ASF_max = max(PopSmetric(index,1));
                            PopSmetric(~index,:) = ASF_max.*ones(sum(~index),1) + sum(PopConScaled(~index,:),2);
                        end
                        %% step 4-3
                        kriging_obj= dacefit(PopDec,PopSmetric,'regpoly0','corrgauss',1*ones(1,D),0.001*ones(1,D),1000*ones(1,D));
                        f_min = min(PopSmetric);
                        %% step 5: Optimization
                        infill_criterion = @(x)Infill_Standard_EI(x, kriging_obj, f_min);
                        best_x=rGA(infill_criterion,Problem);
                        % infill point too close (not used in the paper, 
                        % but this happens sometimes especially for low dimensional problem)
                        if min(sqrt(sum((PopDecT-best_x).^2,2)))<1E-8
                            best_x = rGA(@(x)Infill_Maximal_Distance(x, PopDecT),Problem);
                        end
                        Population = [Population,Problem.Evaluation(best_x)];
                    end   
                end  
            end
        end
    end
end