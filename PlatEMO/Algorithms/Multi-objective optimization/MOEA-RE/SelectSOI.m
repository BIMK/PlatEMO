function [SOI,SOIObjX] = SelectSOI(Population,PopObjX,N,z)
% Select the solutions of interest

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate angles between solutions and ASF
    Angle = acos(1-pdist2(PopObjX,PopObjX,'cosine'));
    Angle(logical(eye(length(Angle)))) = inf;
    ASF = max((PopObjX-z).*sum(PopObjX,2)./PopObjX,[],2);
    
    %% Delete solutions
    Remain = 1 : size(PopObjX,1);
    while length(Remain) > N
        [Dis,x] = min(Angle(Remain,Remain),[],2);
        [~,y]   = min(Dis);
        x = x(y);
        if ASF(Remain(x)) > ASF(Remain(y))
            Remain(x) = [];
        else
            Remain(y) = [];
        end
    end
    SOI     = Population(Remain);
    SOIObjX = PopObjX(Remain,:);
end