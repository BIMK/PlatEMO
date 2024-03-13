function PopDec = EvolALG(Problem,PCheby,Dec,model,IFEs)
% Solution update in ParEGO, where a solution with the best expected
% improvement is re-evaluated

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    Off   = [OperatorGA(Problem,Dec(TournamentSelection(2,size(Dec,1),PCheby),:));OperatorGA(Problem,Dec,{0,0,1,20})];
    N     = size(Off,1);
    EI    = zeros(N,1);
    Gbest = min(PCheby);
    E0    = inf;
    while IFEs > 0
        drawnow('limitrate');
        for i = 1 : N
            [y,~,mse] = predictor(Off(i,:),model);
            s         = sqrt(mse);
            EI(i)     = -(Gbest-y)*normcdf((Gbest-y)/s)-s*normpdf((Gbest-y)/s);
        end
        [~,index] = sort(EI);
        if EI(index(1)) < E0
            Best = Off(index(1),:); 
            E0   = EI(index(1));
        end
        Parent = Off(index(1:ceil(N/2)),:);
        Off    = [OperatorGA(Problem,Parent(TournamentSelection(2,size(Parent,1),EI(index(1:ceil(N/2)))),:));OperatorGA(Problem,Parent,{0,0,1,20})];
        IFEs   = IFEs - size(Off,1);
    end
    PopDec = Best;
end