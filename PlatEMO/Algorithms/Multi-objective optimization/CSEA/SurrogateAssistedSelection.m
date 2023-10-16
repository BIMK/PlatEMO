function Next = SurrogateAssistedSelection(Problem,net,p0,p1,Ref,Input,wmax,tr)
% Surrogate-assisted selection for selecting promising solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    Next  = OperatorGA(Problem,[Input;Ref.decs],{1,15,1,5});
    Label = predict(net,Next);
    a     = tr;
    b     = 1 - tr;
    i     = 0;
    if p0<0.4 || (p1<a&&p0<b)
        while i < wmax
            [~,index] = sort(Label,'descend');
            Input     = Next(index(1:length(Ref)),:);
            Next      = OperatorGA(Problem,[Input;Ref.decs],{1,15,1,5});
            Label = predict(net,Next);
            i = i+size(Next,1);
        end
        Next = Next(Label>0.9,:);
    elseif p0>b && p1<a % Randomly select one to avoid loop
        randindex = randperm(length(Next));
        Next = Next(randindex(1),:);
    elseif p1 > b
        while i<wmax
            [~,index] = sort(Label);
            Input     = Next(index(1:length(Ref)),:);
            Next      = OperatorGA(Problem,[Input;Ref.decs],{1,15,1,5});
            Label = predict(net,Next);
            i = i+size(Next,1);
        end
        Next = Next(Label<0.1,:);
    else
        Next = Next(randi(end),:);
    end
end