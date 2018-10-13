function Next = SurrogateAssistedSelection(Global,net,p0,p1,Ref,Input,wmax,tr)
% Surrogate-assisted selection for selecting promising solutions

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    Next  = Reproduction(Global,Ref,Input);
    Label = net.predict(Next);
    a     = tr;
    b     = 1 - tr;
    i     = 0;
    if p0<0.4 || (p1<a&&p0<b)
        while i < wmax
            [~,index] = sort(Label,'descend');
            Input     = Next(index(1:length(Ref)),:);
            Next      = Reproduction(Global,Ref,Input); 
            Label     = net.predict(Next);
            i = i+size(Next,1);
        end
        Next = Next(Label>0.9,:);
    elseif p0>b && p1<a
        Next = [];
    elseif p1 > b
        while i<wmax
            [~,index] = sort(Label);
            Input     = Next(index(1:length(Ref)),:);
            Next      = Reproduction(Global,Ref,Input); 
            Label     = net.predict(Next);
            i = i+size(Next,1);
        end
        Next = Next(Label<0.1,:);
    else
        Next = [];
    end
end