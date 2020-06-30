function x = LocateWorst(PopObj,W,Region,FrontNo,Z)
% Detect the worst solution in the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Crowd  = hist(Region,1:size(W,1));
    Phi    = find(Crowd==max(Crowd));
    PBI    = CalPBI(PopObj,W,Region,Z,ismember(Region,Phi));
    PBISum = zeros(1,size(W,1));
    for j = 1 : length(PBI)
        PBISum(Region(j)) = PBISum(Region(j)) + PBI(j);
    end
    [~,Phi] = max(PBISum);
    Phih    = find(Region==Phi);
    R       = Phih(FrontNo(Phih)==max(FrontNo(Phih)));
    [~,x]   = max(PBI(R));
    x       = R(x);
end