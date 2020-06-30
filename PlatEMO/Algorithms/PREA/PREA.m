function PREA(Global)
% <algorithm> <P>
% Promising-region based EMO algorithm

%------------------------------- Reference --------------------------------
% J. Yuan, H. Liu, F. Gu, Q. Zhang, and Z. He, Investigating the properties
% of indicators and an evolutionary many-objective algorithm based on a
% promising region, IEEE Transactions on Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

%% Generate the random population
Population   = Global.Initialization();
Zmin         = min(Population.objs,[],1);

%% shift the objective space to R+
PopObj               = Population.objs;
PopObj               = PopObj - repmat(Zmin,Global.N,1) + 1e-6;

%% calculate the ratio based indicator matrix
IMatrix              = ones(Global.N,Global.N); 
for i=1:1:Global.N
    Fi               = PopObj(i,:);
    % calculate ratio based indicator value of each individual
    Ir               = PopObj./repmat(Fi,Global.N,1) - 1;  
    InvertIr         = repmat(Fi,Global.N,1)./PopObj - 1;
    MaxIr            = max(Ir,[],2);
    MinIr            = max(InvertIr,[],2);
    DomInds          = find(MaxIr<=0);
    MaxIr(DomInds)   = -MinIr(DomInds);
    IMatrix(i,:)     = MaxIr';
    IMatrix(i,i)     = Inf;
end

%% Optimization
while Global.NotTermination(Population)
    MatingPool       = MatingStrategy(IMatrix);
    Offspring        = GAhalf(Population(MatingPool));
    Zmin             = min([Zmin;Offspring.objs],[],1);
    [Population,IMatrix] = PREA_Update([Population,Offspring],Global.N,Zmin);
end
end