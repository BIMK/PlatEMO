function [grade,X,SR,improve,AppSet] = LocalSearch3(Global,X,SR,improve,AppSet)
% Local Search 3 of MTS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    grade   = 0;
    SearchL = Global.lower;
    SearchU = Global.upper;
    Disp    = (SearchU-SearchL)/10;
    Best    = X;
    while any(Disp<1e-2)
        for i = randperm(length(SR))
            value = SearchL(i) : Disp(i) : SearchU(i);
            Decs  = repmat(Best.dec,length(value),1);
            Decs(:,i) = value;
            Y     = INDIVIDUAL(Decs);
            for j = 1 : length(Y)
                % 'grade', 'SR' and 'improve' will not be updated in local
                % search 3
                [~,~,AppSet] = Grading(Y(j),Best,grade,improve,AppSet,Global.N);
                if all(Y(j).obj<=Best.obj)
                    Best = Y(j);
                end
            end
            SearchL(i) = max(Best.dec(i)-2*Disp(i),Global.lower);
            SearchU(i) = min(Best.dec(i)+2*Disp(i),Global.upper);
            Disp(i)    = (SearchU(i)-SearchL(i))/10;
        end
    end
    X = Best;
end