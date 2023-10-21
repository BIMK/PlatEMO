function SDE = CalSDE(PopObj)
% Calculate the value of SDE of each solution

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N      = size(PopObj,1);   
    Zmin   = min(PopObj,[],1);
    Zmax   = max(PopObj);
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(Zmax,N,1)-repmat(Zmin,N,1));
    SDE    = zeros(1,N);
    for i = 1 : N
        SPopuObj = PopObj;
        Temp     = repmat(PopObj(i,:),N,1);
        Shifted  = PopObj < Temp;
        SPopuObj(Shifted) = Temp(Shifted);                                    
        Distance  = pdist2(real(PopObj(i,:)),real(SPopuObj));
        [~,index] = sort(Distance,2);
        Dk = Distance(index(floor(sqrt(N))+1)); % Dk denotes the distance of solution i and its floor(sqrt(N)+1)-th nearest neighbour
        SDE(i)=1./(Dk+2);
    end
end