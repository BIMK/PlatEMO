function [grade,improve,AppSet] = Grading(X,oldX,grade,improve,AppSet,popsize)
% Check if the new solution can be added to the approximation set

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Decide whether to put the new solution into approximation set
    if isempty(AppSet)
        AppSet = X;
        flag   = false;
    else
        compare = AppSet.objs - repmat(X.obj,length(AppSet),1);
        if any(all(compare<=0,2))
            flag = false;
        else
            dominated = all(compare>=0,2);
            AppSet(dominated) = [];
            AppSet = [AppSet,X];
            flag   = true;
        end
    end
    % If the archive contains too many solutions, delete some. The original
    % algorithm does not limit the size of archive, so that the size of
    % archive will increase unrestrainedly
    if length(AppSet) > 10*popsize
        AppSet = AdjustAppSet(AppSet,5*popsize);
    end

    %% Update the information about the solution
    if flag
        grade = grade + 9;
    end
    if sum(X.obj<oldX.obj) > sum(X.obj>oldX.obj)
        grade   = grade + 2;
        improve = true;
    end
end