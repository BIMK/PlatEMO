function ArcPop = UpdateArc(Population,N)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [FrontNo,MaxFNo] = NDSort([Population.objs,sum(max(0,Population.cons),2)],1);
    Temp1 = FrontNo==1;

    Population = Population(Temp1==1);

    Temp = find(sum(max(0,Population.cons),2)>0);
    Population = Population(Temp);

    if length(Population)<N 
        ArcPop = Population;
    else
        Zmax = max(Population.objs,[],1);
        Next(1:size(Population,2)) = true;
        % Select the solutions in the last front
        Delete = LastSelection(Population(Next).objs,-sum(max(0,Population.cons),2),sum(Next)-N,Zmax);
        Temp = find(Next);
        Next(Temp(Delete)) = false;
        ArcPop = Population(Next);
    end
end


function Delete = LastSelection(PopObj,PopCons,K,Zmax)
% Select part of the solutions in the last front

    [N,M]  = size(PopObj);
    PopObj = (PopObj-repmat(Zmax,N,1))./(repmat(min(PopObj),N,1)-repmat(Zmax,N,1)- 1e-10);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine = 1 - pdist2(PopObj,PopObj,'cosine');
    Cosine = Cosine.*(1-eye(size(PopObj,1)));

    %% Environmental selection
    Delete = false(1,N);
    % Select K solutions one by one
    while sum(Delete) < K
        [Jmin_row,Jmin_column] = find(Cosine==max(max(Cosine)));
        j = randi(length(Jmin_row));
        Temp_1 = Jmin_row(j);
        Temp_2 = Jmin_column(j);
        if (PopCons(Temp_1)<PopCons(Temp_2)) ||(PopCons(Temp_1)==PopCons(Temp_2) && rand<0.5)
            Delete(Temp_1) = true;
            Cosine(:,Temp_1)=0;
            Cosine(Temp_1,:)=0;
        else
            Delete(Temp_2) = true;
            Cosine(:,Temp_2)=0;
            Cosine(Temp_2,:)=0;
        end
    end
end