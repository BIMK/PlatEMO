function Next = SurrogateAssistedSelectionPC(Problem,net,error1,error2,Input,wmax,Pa,D,flag,delta)
% Surrogate-assisted selection for selecting promising solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Pa    = Pa(:,1:D);
    Next  = OperatorGA(Problem,[Input;Pa(:,1:D)],{1,15,1,5});
    Label = net.lastpredict(Next,D,Pa,flag);
    lnum  = size(Pa,1);
    i     = 1;
    wmax  = floor(wmax/lnum);
    GoodNext  = zeros(wmax,D);
    GoodLabel = zeros(wmax,1);
    if error1 < 1-delta
        while i <= wmax
            [~,index] = sort(Label,'descend');
            GoodNext(i,:)  = Next(index(1),:);
            GoodLabel(i,:) = Label(index(1));
            Input  = Next(index(lnum),:);
            Parent = [Input;Pa];
            Parent = Parent(randperm(size(Parent,1)),:);
            Next   = OperatorGA(Problem,Parent,{1,15,1,5});
            Label  = net.lastpredict(Next,D,Pa,flag);
            i = i+1;
        end
        if (sum(GoodLabel >= 0.95) == 0) || (sum(GoodLabel >= 0.95) > floor(lnum/2))
            [~,index] = sort(GoodLabel(:,end),'descend');
            Next = GoodNext(index(1:floor(lnum/2)),1:D);
        else
            Next = GoodNext(GoodLabel >= 0.95,:);
        end

    elseif error2 < 1-delta
        while i <= wmax
            [~,index] = sort(Label);
            GoodNext(i,:)  = Next(index(1),:);
            GoodLabel(i,:) = Label(index(1));
            Input  = Next(index(lnum),:);
            Parent = [Input;Pa];
            Parent = Parent(randperm(size(Parent,1)),:);
            Next   = OperatorGA(Problem,Parent,{1,15,1,5});
            Label  = net.lastpredict(Next,D,Pa,flag);
            i = i+1;
        end
        if (sum(GoodLabel <= -0.95) == 0) || (sum(GoodLabel <= -0.95) > floor(lnum/2))
            [~,index] = sort(GoodLabel(:,end));
            Next = GoodNext(index(1:floor(lnum/2)),1:D);
        else
            Next = GoodNext(GoodLabel <= -0.95,:);
        end
    else
        Next = Next(randi(end),:);
    end
end