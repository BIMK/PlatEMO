function Population = PUpdate(Population,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao.zhang.cn@gmail.com)

    %% Data Preprocess
    global V
    PopObj = Population.objs ; PopCon = Population.cons;
    zmin   = min(PopObj,[],1); zmax = max(PopObj,[],1);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    NV     = size(V,1);
    NP     = size(PopObj,1);

    %% Associate each solution to a reference vector
    Angle  = acos(1-pdist2(PopObj,V,'cosine'));
    Pindex = true(1,NP);
    Vindex = true(1,NV);
    
    %% Select one solution for each reference vector
    while any(Vindex)
        [~,associate] = min(Angle(Pindex,Vindex),[],2);
        Pexist = find(Pindex==1); 
        Vexist = find(Vindex==1); 
        for i = unique(associate)'
            current = find(associate==i);
            if isscalar(current)
                best = 1;
            else
                best = scalar_select(PopObj(Pexist(current),:),PopCon(Pexist(current),:),V(Vexist(i),:));
            end
            Pindex(Pexist(current(best))) = 0;
            Vindex(Vexist(i)) = 0;
        end
    end

    %% Select remain solutions
    if N > NV
        Pexist = find(Pindex==1);
        index  = randperm(length(Pexist),N-NV);
        Pindex(Pexist(index)) = 0;
    end
    Next = Pindex==0;
    
    %% Population for next generation
    Population = Population(Next);
end

function best = scalar_select(PopObj,PopCon,lamda)
    g    = max(PopObj.*lamda,[],2);
    CV   = sum(max(0,PopCon),2);
    best = find(CV==min(CV,[],1));
    if length(best) > 1
        [~,best_] = min(g(best));
        best = best(best_);
    end
end