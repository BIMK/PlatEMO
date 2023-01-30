function [OffDec,OffMask] = Operator_pvfv(Problem,ParentDec,ParentMask,pv,fv,delta)
% The operator of MSKEA guided by pv and sv

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lei Chen

    %% Parameter setting
    [N,D]       = size(ParentDec);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);
    
    %% Crossover for mask
    OffMask = Parent1Mask;
    for i = 1 : N/2
        if rand < 0.5
            index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
            index = index(TS(-pv(index)));
            OffMask(i,index) = 0;
        else
            index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
            index = index(TS(pv(index)));
            OffMask(i,index) = Parent2Mask(i,index);
        end
    end
    
    %% Mutation for mask
    if rand<(1-delta)
        f_vector=fv;
        f_vector(fv>0)=1;
        for i=1:N/2
            index=find(OffMask(i,:)~=f_vector);
            if rand<0.5
                index=index(TS(-fv(index)));
                OffMask(i,index)=1;
            else
                index=index(TS(fv(index)));
                OffMask(i,index)=0;
            end
        end
    else
        for i = 1 : N/2
            if rand < 0.5
                index = find(OffMask(i,:));
                index = index(TS(-pv(index)));
                OffMask(i,index) = 0;
            else
                index = find(~OffMask(i,:));
                index = index(TS(pv(index)));
                OffMask(i,index) = 1;
            end
        end
    end
    
    %% Crossover and mutation for dec
    if any(Problem.encoding~=4)
        OffDec = OperatorGAhalf(Problem,ParentDec);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(N/2,D);
    end
end

function index = TS(pv)
% Binary tournament selection

    if isempty(pv)
        index = [];
    else
        index = TournamentSelection(2,1,pv);
    end
end