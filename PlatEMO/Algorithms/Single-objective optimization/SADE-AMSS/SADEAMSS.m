classdef SADEAMSS < ALGORITHM
% <2023> <single> <real/integer> <expensive>
% Surrogate-assisted differential evolution with adaptive multi-subspace search
% K    ---  20 --- Maximum number of subspaces in a generation
% maxd --- 100 --- Maximum number of variables in a subspace
% Gm   ---   5 --- Maximum iterations of subspace optimization

%------------------------------- Reference --------------------------------
% H. Gu, H. Wang, and Y. Jin. Surrogate-assisted differential evolution
% with adaptive multi-subspace search for large-scale expensive
% optimization. IEEE Transcations on Evolutionary Computation, 2023, 27(6):
% 1765-1779.
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
            %% Parameter setting
            [K,maxd,Gm] = Algorithm.ParameterSet(20,100,5);

            %% Initialization
            % Initialization of DSD
            fsum  = 0;
            t     = 1;
            Emean = 1;
            DSD   = 100000;
            tstop = 0;
            switchflag = 0;
            tes   = 50;
            tr    = 500;
            beta  = 2;
            Fes   = 0;
            % Initialization of population
            Ns        = 200;
            Archive   = Problem.Evaluation(repmat(Problem.upper-Problem.lower,Ns,1).*UniformPoint(Ns,Problem.D,'Latin')+repmat(Problem.lower,Ns,1));
            bestvalue = min(Archive.objs);
            popt      = repmat(Problem.upper-Problem.lower,Problem.N,1).*UniformPoint(Problem.N,Problem.D,'Latin') + repmat(Problem.lower,Problem.N,1);
            % Lower and upper boundary points
            slu      = [Problem.lower;Problem.upper];    
            pop_mean = 0;
            vec_pop  = eye(Problem.D);
            strflag  = 2;

            %% Optimize
            while Algorithm.NotTerminated(Archive)
                if mod(Fes+1,strflag) == 0  % Mapping subspace
                    [~,pop_mean] = zscore(popt);
                    popt = popt-repmat(pop_mean,size(popt,1),1);
                    covx = cov(popt);
                    [vec_pop,lamd] = pcacov(covx);  % Descending order
                    a    = lamd(lamd > 1);
                    d_1  = size(a,1);
                    popt = popt * vec_pop;
                    slu  = [Problem.lower;Problem.upper];
                    slu  = slu - repmat(pop_mean,2,1);
                    slu  = sort(slu*vec_pop,1);
                end
                if mod(Fes+1,strflag) ~= 0  % Original subspace
                    slu      = [Problem.lower;Problem.upper];
                    pop_mean = 0;
                    vec_pop  = eye(Problem.D);
                end
                k     = 1;
                [~,b] = min(Archive.objs);  % bestX for DE
                bestX = Archive(b).dec;
                bestX = bestX - repmat(pop_mean,size(bestX,1),1);
                bestX = bestX*vec_pop;
                while k <= K
                    d = randi(min(Problem.D,maxd));  % Number of decision variables in the k-th subspace
                    if mod(Fes+1,strflag) == 0
                        if d > d_1         % Case_1
                            q2   = d - d_1;
                            col_rand2 = randperm(Problem.D-d_1);
                            col2 = col_rand2(1:q2);
                            col2 = d_1 + col2(1:q2);
                            col  = [1:d_1,col2];
                            randIndex = randperm(size(col,2));
                            col = col(:,randIndex);
                        elseif d <= d_1   % Case_2
                            col_rand1 = randperm(d);
                            col = col_rand1(1:d);
                        end
                    end
                    if mod(Fes+1,strflag) ~= 0
                        col_rand = randperm(Problem.D);
                        col = col_rand(1:d);
                    end
                    popk = popt(:,col);
                    tsn  = 2*d;     % Size of the training set
                    Xtrain_rand    = randperm(length(Archive));
                    ip_Xtrain_rand = Xtrain_rand(1:tsn);
                    x_kth_trains   = Archive(ip_Xtrain_rand).decs;
                    x_kth_trains   = x_kth_trains - repmat(pop_mean,size(x_kth_trains,1),1);
                    x_kth_trains   = x_kth_trains * vec_pop;
                    x_kth_train    = x_kth_trains(:,col);
                    arc_train      = Archive(ip_Xtrain_rand).objs;
                    k_train_point{k}     = x_kth_train;
                    [lambda{k},gamma{k}] = RBF(x_kth_train,arc_train,'cubic');  % RBFN model
                    % Differential evolution
                    XG  = popk;
                    XGf = RBF_eval(XG,x_kth_train,lambda{k},gamma{k},'cubic');
                    g   = 1;
                    while g <= Gm
                        XG_next  = Operator(slu(:,col),XG,bestX(1,col));
                        f_next_G = RBF_eval(XG_next,x_kth_train,lambda{k},gamma{k},'cubic');
                        now_popf = [XGf;f_next_G];
                        [~,y]    = sort(now_popf);
                        now_pop  = [XG;XG_next];
                        % Select best Np to the next g
                        XG  = now_pop(y(1:Problem.N),:);
                        XGf = now_popf(y(1:Problem.N));
                        g = g + 1;
                    end
                    spk   = XGf(1);     % Best result of k-th subspace
                    xbspk = XG(1,:);
                    popt(:,col)  = XG;  % Update the population
                    colfff{k}    = col;
                    Dspkfff(k)   = d;
                    bestspk(k)   = spk;
                    bestxbspk{k} = xbspk;
                    k = k + 1;
                end
                xbest_t = bestX;
                [value_minxbest,pos_minxbest] = min(bestspk);
                kthbestx = bestxbspk{pos_minxbest};
                bestcol  = colfff{pos_minxbest};
                popt_kth = popt(1,bestcol);
                out_popt_kth = RBF_eval(popt_kth,k_train_point{pos_minxbest},lambda{pos_minxbest},gamma{pos_minxbest}, 'cubic');
                if out_popt_kth <= value_minxbest
                    xbest_t(:,bestcol) = popt_kth;
                else
                    xbest_t(:,bestcol) = kthbestx;
                end
                Fes = Fes + 1;
                % Date space transformation
                xpz = xbest_t/vec_pop;
                xp  = xpz + repmat(pop_mean,size(xpz,1),1);
                xp  = min(max(xp,Problem.lower),Problem.upper);
                if mod(Fes,strflag) == 0
                    popt = popt/vec_pop;    % Date space transformation
                    popt = popt + repmat(pop_mean,size(popt,1),1);
                    popt = min(max(popt,Problem.lower),Problem.upper);
                end
                Archive   = [Archive,Problem.Evaluation(xp)];
                bestvalue = min(min(Archive.objs),bestvalue);
                % Adaptive switching strategy
                if Fes > 2
                    bestvalue1 = log10(min(Archive(1:Fes+Ns-2).objs));
                    bestvalue2 = log10(min(Archive(1:Fes+Ns-1).objs));
                    bestvalue3 = log10(min(Archive(1:Fes+Ns).objs));
                    change1    = bestvalue2 - bestvalue1;
                    change2    = bestvalue3 - bestvalue2;
                    f_sd(t)    = (change2-change1)/2;
                    t = t + 1;
                    if t == 51
                        for i = 1 : tes
                            fsum = fsum + abs(f_sd(i));
                        end
                        Emean = fsum/tes;
                    end
                    f_sd_mn = f_sd/Emean;
                    if t > tr+tstop
                        DSD = sum(abs(f_sd_mn(t-tr:t-1)));
                    end
                    if DSD < beta
                        tstop      = Fes + tr;
                        switchflag = switchflag + 1;
                    end
                    if switchflag == 1
                        strflag = 4;
                        DSD     = 100000;
                    end
                    if switchflag == 2
                        strflag = 1000000000;
                        DSD     = 100000;
                    end
                end
            end
        end
    end
end