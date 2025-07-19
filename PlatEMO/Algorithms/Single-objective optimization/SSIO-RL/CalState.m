function State = CalState(Pdec, Pobj, Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------


    N           = length(Pobj);    
    [Pobj,rank] = sort(Pobj);
    Pdec        = Pdec(rank,:);
    distance    = pdist2(Pdec,Pdec(1,:));
    ccc         = cov(distance,Pobj);
    state1      = ccc(1,2)/(sqrt(ccc(1,1))*sqrt(ccc(2,2)));
    state2      = 1 - (max(Pobj)-min(Pobj))/max(Pobj);
    state3      = 1 - (max(Pobj)-mean(Pobj))/max(Pobj);
    state4      = std(Pobj)/std([repmat(max(Pobj),N/2,1);repmat(min(Pobj),N/2,1)]);
    state5      = mean(pdist2(Pdec(1,:),Pdec))/pdist2(max(Pdec),min(Pdec));
    state6      =  1 - Problem.FE/Problem.maxFE;
    State       = [state1;state2;state3;state4;state5;state6];
end

