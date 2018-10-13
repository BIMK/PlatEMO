function [zmin,zmax,Record,Mark] = UpdateReferencePoints(PopObj,zmin,zmax,Record,Mark,alpha,epsilon)
% Update the minimum and maximum values of each objective for normalization

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    z      = min(PopObj,[],1);
    znad   = max(PopObj,[],1);
    zmin   = min(zmin,z);
    Record = [Record(2:end,:);znad];
    v      = Record(end-1,:) - znad;
    mark   = false(1,length(zmax));
    if max(v) > alpha
        zmax = znad;
    else
        for i = 1 : length(zmax)
            if abs(zmax(i)-zmin(i)) < epsilon
                zmax(i) = max(zmax);
                mark(i) = true;
            elseif znad(i) > zmax(i)
                zmax(i) = 2*znad(i) - zmax(i);
                mark(i) = true;
            elseif v(i)==0 && ~any(Mark(:,i))
                zmax(i) = (zmax(i)+max(Record(:,i)))/2;
                mark(i) = true;
            end
        end
    end
    Mark = [Mark(2:end,:);mark];
end