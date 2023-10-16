function Select = subsetSelection(Obj,Objhat,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    Len = length(Obj);
    Obj = [Obj;Objhat];
    
    %% Select the representative objective vector set
    LpNormD = pdist2(Obj,Obj);
    Select = false(1,size(Obj,1));
    Select(1:Len) = true;
    % Greedy inclusion distance-based subset slection
    while sum(Select) < N
        Remain   = find(~Select);
        [~, rho] = max(min(LpNormD(Remain,Select),[],2));
        Select(Remain(rho)) = true;
    end
    Select = Select(Len+1:end);
end