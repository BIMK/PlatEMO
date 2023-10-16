classdef MCCMO< ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Multi-population coevolutionary constrained multi-objective optimization

%------------------------------- Reference --------------------------------
% J. Zou, R. Sun, Y. Liu, Y. Hu, S. Yang, J. Zheng, and K. Li, A
% multi-population evolutionary algorithm using new cooperative mechanism
% for solving multi-objective problems with multi-constraint, IEEE
% Transactions on Evolutionary Computation, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            Pop0 = Problem.Initialization();
            Fit0 = CalFitness(Pop0.objs);
            totalcon = size(Pop0(1,1).con,2);
            Obj0(1) = sum(sum(Pop0.objs,1));
            Pops = [];
            gen = 2;
            change_threshold = 1e-2;
            constant = 1;
            color = [0,0,1;1,0,0;0,1,0;1,0,1;1,1,0;0,1,1;1,1,1;0.5,0.3,0.2;0.7,0.8,0.1;0.4,0.3,0.9;1,0.1,0.5;0.2,0.8,1;0.55,0.28,0.81;0.1,0.8,0.3;0.5,0.3,0.7];
            dw  = 0;
            Off = Pop0;
            minep = 1e-5;
            index = 1;
            for i = 1:totalcon
                Pops{1,i} = Problem.Initialization();
                Pops{2,i} = CalFitness(Pops{1,i}.objs);
                Pops{4,i} = i;
                processedcon{index} = i;
                index = index+1;
                Pops{3,i}(1) = sum(sum(Pops{1,i}.objs,1));
                Pops{5,i} = 0;
                if dw
                    %                         Pops{6,i} = color(i,:);
                    Pops{6,i} = [rand(1),rand(1),rand(1)];
                end
                Pops{7,i} = 1;
                Pops{8,i} = i;
                Off = [Off,Pops{1,i}];
            end
            processedcon{index} = totalcon+1;
            index = index+1;
            
            UPF = 0;
            combine = 0;
            reinit = 0;
            arch = [];
            [arch,Fit_arch] = EnvironmentalSelection([arch,Off],Problem.N,totalcon+1,totalcon);
            numberpop(1) = 2;
            evad(1) = totalcon*100+100;
            ctime = 0;
            flag777=0;
            act=0;
            picked = [];
            Feas = [];
            IGDV = [];
            dUPF = 0;
                        
            %% Optimization
            while Algorithm.NotTerminated(arch)
                Off = [];
                Offtemp = [];
                
                if contains( class(Problem),'RWMOP') || contains(class(Problem),'LIRCMOP') || contains(class(Problem),'DOC') || contains(class(Problem),'DASCMOP')
                    for i = 1:size(Pops,2)
                        Offtemp = [];
                        if size(Pops,2) > 1 && ~Pops{7,i}
                            MatingPool = TournamentSelection(2,Problem.N,Pops{2,i});
                            Offtemp  = OperatorDE(Problem,Pops{1,i}(1:Problem.N/2),Pops{1,i}(MatingPool(1:end/2)),Pops{1,i}(MatingPool(end/2+1:end)));
                        end
                        Off = [Off,Offtemp];
                    end
                    if ~UPF || ~dUPF
                        MatingPool2 = TournamentSelection(2,Problem.N,Fit0);
                        Offtemp = OperatorDE(Problem,Pop0(1:end/2),Pop0(MatingPool2(1:end/2)),Pop0(MatingPool2(end/2+1:end)));
                        Off = [Off,Offtemp];
                    elseif dUPF && reinit
                        MatingPool1 = TournamentSelection(2,Problem.N/2,Fit_arch);
                        MatingPool2 = TournamentSelection(2,Problem.N/2,Fit0);
                        Offtemp = OperatorDE(Problem,Pop0(1:end/2),Pop0(MatingPool2),arch(MatingPool1));
                        Off = [Off,Offtemp];
                    end
                    MatingPool1 = TournamentSelection(2,2*Problem.N,Fit_arch);
                    Offtemp = OperatorDE(Problem,arch,arch(MatingPool1(1:end/2)),arch(MatingPool1(end/2+1:end)));
                    Off = [Off,Offtemp];
                else
                    for i = 1:size(Pops,2)
                        MatingPool = TournamentSelection(2,floor(Problem.N),Pops{2,i});
                        if size(Pops,2) > 1 && ~Pops{7,i}
                            Offtemp  = OperatorGAhalf(Problem,Pops{1,i}(MatingPool));
                        end
                        Off = [Off,Offtemp];
                    end
                    if ~UPF || ~dUPF
                        MatingPool2 = TournamentSelection(2,Problem.N,Fit0);
                        Offtemp  = OperatorGAhalf(Problem,Pop0(MatingPool2));
                        Off = [Off,Offtemp];
                    elseif dUPF && reinit
                        MatingPool1 = TournamentSelection(2,Problem.N/2,Fit_arch);
                        MatingPool2 = TournamentSelection(2,Problem.N/2,Fit0);
                        Offtemp  = OperatorGAhalf(Problem,[Pop0(MatingPool2),arch(MatingPool1)]);
                        Off = [Off,Offtemp];
                    end
                    MatingPool1 = TournamentSelection(2,Problem.N,Fit_arch);
                    Offtemp  = OperatorGA(Problem,arch(MatingPool1));
                    Off = [Off,Offtemp];
                end

                [arch,Fit_arch,selected] = EnvironmentalSelection([Off,arch],Problem.N,totalcon+1,totalcon);
                [Pop0,Fit0] = EnvironmentalSelection([Pop0,Off],Problem.N,0,totalcon);
                Obj0(gen) = sum(sum(abs(Pop0.objs),1));
                if (dw && reinit) || (dw && ~dUPF)
                    Draw(Pop0.objs,'sk','Markeredgecolor',[.5 .5 .5],'Markerfacecolor',[.8 .5 .5]);
                end
                if size(Pops,2) > 1
                    for i = 1:size(Pops,2)
                        if dw
                            Draw(Pops{1,i}.objs,'o','Markeredgecolor','black','Markerfacecolor',Pops{6,i});
                        end
                        [Pops{1,i},Pops{2,i}] = EnvironmentalSelection([Pops{1,i},Off],Problem.N,Pops{4,i},totalcon,0);
                        Pops{3,i}(gen) = sum(sum(Pops{1,i}.objs,1));
                        
                        if size(Pops,2) > 1
                            Pops{5,i} = is_stable(Pops{3,i},gen,Pops{1,i},Problem.N,change_threshold,constant,Problem.M,~Pops{7,i});
                        else
                            Pops{5,i} = 1;
                        end
                        if Pops{5,i} && Pops{7,i}
                            Pops{7,i} = 0;
                            Pops{5,i} = 0;
                        end
                    end
                end
                
                if UPF == 0
                    UPF = is_stable(Obj0,gen,Pop0,Problem.N,change_threshold,constant,Problem.M,0);
                else
                    dUPF = is_stable(Obj0,gen,Pop0,Problem.N,change_threshold,constant,Problem.M,0);
                end
                if dUPF
                    pop_cons = Pop0.cons;
                    cv = overall_cv(pop_cons);
                    FR = size(find(cv <= 0),1) / Problem.N;
                    if FR
                        reinit = 1;
                    else
                        reinit = 0;
                    end
                end

                if size(Pops,2) > 1
                    i = 1;
                    AllPops0 = [];
                    for j = 1:size(Pops,2)
                        AllPops0 = [AllPops0,Pops{1,j}];
                    end
                    if UPF && dUPF
                        AllPops0 = [AllPops0,arch];
                    end
                    [FrontNo,~] = NDSort(AllPops0.objs,inf);
                    if UPF && dUPF
                        Minindex = min(FrontNo(size(Pops,2)*Problem.N+1:end));
                        for j = 1:size(Pops,2)
                            if max(FrontNo((j-1)*Problem.N+1:j*Problem.N))< Minindex
                                Pops{5,j} = 1;
                            end
                        end
                    end
                    while i <= size(Pops,2)
                        if Pops{5,i}
                            Minindex = min(FrontNo((i-1)*Problem.N+1:i*Problem.N));
                            for j = 1:size(Pops,2)
                                if i ~= j
                                    if max(FrontNo((j-1)*Problem.N+1:j*Problem.N))< Minindex
                                        Pops{5,j} = 1;
                                    end
                                end
                            end
                        end
                        i = i+1;
                    end
                    
                    merge_index = [];
                    merge_pop = [];
                    merge_con = [];
                    domine = 0;
                    for p = 1:size(Pops,2)
                        if Pops{5,p} == 1
                            conp =  Pops{1,p}.cons;
                            conp = max(0,conp(:,Pops{4,p}));
                            if size(find(conp <= 0),1)
                                merge_pop = [merge_pop,Pops{1,p}];
                                merge_index = [merge_index,p];
                                merge_con = [merge_con,Pops{4,p}];
                            else
                                Pops{1,p} = Problem.Initialization();
                                Pops{2,p} = CalFitness(Pops{1,p}.objs,Pops{1,p}.cons,Pops{4,p});
                                Pops{9,p} = 0;
                            end
                        end
                    end
                    if size(merge_index,2) > 1
                        selected = merge_index(1);
                        Pops{4,selected} = merge_con;
                        [Pops{1,selected},Pops{2,selected}] = EnvironmentalSelection([AllPops0,arch],Problem.N,Pops{4,selected},totalcon);
                        Pops{5,selected} = 0;
                        Pops{3,selected}(gen) = sum(sum(Pops{1,selected}.objs,1));
                        Pops{7,selected} = 0;
                        Pops{8,selected} = index;
                        processedcon{index} = merge_con;
                        index = index+1;
                        k = size(Pops,2);
                        while k
                            if Pops{5,k} == 1
                                Pops(:,k) = [];
                            end
                            k = k-1;
                        end
                        
                    end
                    
                    if size(Pops,2) == 1
                        [arch,Fit_arch] = EnvironmentalSelection([arch,Pops{1,1}],Problem.N,totalcon+1,totalcon);
                    end
                end
                gen = gen+1;
                if Problem.FE >= Problem.maxFE
                    a =1 ;
                end
            end
        end
    end
end

function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

function result = is_stable(Objvalues,gen,Population,N,change_threshold,constant,M,isND)
    result = 0;
    if isND
        [FrontNo,~]=NDSort(Population.objs,size(Population.objs,1));
        NC=size(find(FrontNo==1),2);
    else
        NC =N;
    end
    max_change = abs(Objvalues(gen)-Objvalues(gen-1));

    if NC == N
        chth = change_threshold * abs((Objvalues(gen) )/( N* M))*10^(M-2);
        if max_change <= chth
            result = 1;
        end
    end
end