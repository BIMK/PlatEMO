function [Population,Archive] = EnvironmentalSelection(Population,Offspring,Archive,MaxSize,z,znad,theta,CSAObjs,epsilon)
% Environmental selection of CMaDPPs

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    PopOff  = [Offspring,Population,Archive];
    PopObjs = PopOff.objs;
    PopCons = PopOff.cons;
    
    ND = NDSort(PopObjs,PopCons,inf);
    Population1 = PopOff(ND==1);
    
    IMatrix = Indicator_Cal(PopOff,z);
    IrFitness1  = min(IMatrix,[],2);
    Norm_IrFitness1 = IrFitness1;

    Population2 = PopOff(Norm_IrFitness1 >= 0);
    
    Population = [Population1,Population2];
    
    CAObj=Population.objs;
    [~,ia,~] = unique(CAObj,'rows');
    Population = Population(ia);
    N = length(Population);
    
    if N <= MaxSize
        return;
    end
    
    IMatrix = Indicator_Cal(Population,z);
    IrFitness  = min(IMatrix,[],2);
    Norm_IrFitness = IrFitness;
    
    L = Norm_IrFitness*Norm_IrFitness';
    
    CAObj = Population.objs;
    CAObj2 = (CAObj-repmat(z,N,1))./(repmat(znad,N,1)-repmat(z,N,1));
   
    D = pdist2(CAObj2,CAObj2,'cosine');

    D(1:size(D,1)+1:end) = 0;
    
    H=(1-theta)*exp(-D);

    KM = H.*L;
    
    L = decompose_kernel(KM);
    Choose = sample_dpp(L,MaxSize);
    Population = Population(Choose);
    
    PopCons = PopOff.cons;
    CV = sum(max(0,PopCons),2);
    FeasibleInd = CV <= epsilon;
    Archive = UpdateArchive(PopOff(FeasibleInd),MaxSize,z,znad,CSAObjs,theta,epsilon);
end