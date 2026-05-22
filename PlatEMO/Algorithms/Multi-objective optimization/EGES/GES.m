function [PopDec,PopObj] = GES(A1,Model,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    typeOfGroups   = 2;
    numberOfGroups = 4;

    [V, Problem.N] = UniformPoint(Problem.N,Problem.M); % Reference Points
    
    Zmin = min(A1.objs,[],1);
    if size(A1.decs,1) >= Problem.N
        Next = NSGAIIIEnvironmentalSelection(A1,Problem.N,V,Zmin); 
    end
    
    Population = Next;
    Zmin       = min(Population(all(Population.cons<=0,2)).objs,[],1); % Ideal point
    
    w    = 1;
    wmax = 20;
    while w <= wmax
        MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));                        % Tournament Selection
        Offspring  = GESoperating(Problem, Population(MatingPool), numberOfGroups, typeOfGroups, Model, A1);    % GES
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);                             % Ideal point
        Population = NSGAIIIEnvironmentalSelection([Population,Offspring],Problem.N,V,Zmin);                % GLMO_NSGAIIIEnvironmentalSelection
        w          = w+1;
    end
    PopDec = Population.decs;
    PopObj = Population.objs;
end