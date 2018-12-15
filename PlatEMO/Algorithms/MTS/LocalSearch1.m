function [grade,X,SR,improve,AppSet] = LocalSearch1(Global,X,SR,improve,AppSet)
% Local Search 1 of MTS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if ~improve
        SR = SR / 2;
        if all(SR<1e-8)
            SR = (Global.upper-Global.lower).*(rand(1,length(SR))/10+0.4);
        end
    end
    improve = false;
    grade   = 0;
    for i = randperm(length(SR))
        old_X  = X;
        dec    = X.dec;
        dec(i) = dec(i) + SR(i)*(rand*2-1);
        X      = INDIVIDUAL(dec);
        [grade,improve,AppSet] = Grading(X,old_X,grade,improve,AppSet,Global.N);
        if all(old_X.obj<=X.obj)
            dec    = old_X.dec;
            dec(i) = dec(i) - 0.5*SR(i)*(rand*2-1);
            X      = INDIVIDUAL(dec);
            [grade,improve,AppSet] = Grading(X,old_X,grade,improve,AppSet,Global.N);
            if all(old_X.obj<=X.obj)
                X = old_X;
            end
        end
    end
end