function [OffDec,OffMask] = Operator_sv(Problem,ParentDec,ParentMask,sv)
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
    rate0   = sv;           % The probability that  0 inverts to  1
    rate1   = 1-rate0;      % The probability that  1 inverts to  0

    for i = 1 : N/2
        diff = find(Parent1Mask(i,:)~=Parent2Mask(i,:));
        temp_rate1=rate1(diff);
        temp_rate0=rate0(diff);
        rate = zeros(1,length(diff));
        rate(logical(OffMask(i,diff)))  = temp_rate1(logical(OffMask(i,diff)));
        rate(logical(~OffMask(i,diff))) = temp_rate0(logical(~OffMask(i,diff)));
        exchange  = rand(1,length(diff)) < rate;
        OffMask(i,diff(exchange))=~OffMask(i,diff(exchange)); 
    end

    %% Mutation for mask
    Mutation_p=1/D;                     % Probability of mutation
    Mu_exchange=rand(N/2,D)<Mutation_p; %The decision variables less than 1/D will be mutated
    for i=1:N/2
        if sum(Mu_exchange(i,:))
            subscript =find(Mu_exchange(i,:)==1);
            rate=zeros(1,size(subscript,2));
            rate(logical(OffMask(i,subscript)))=rate1(subscript(logical(OffMask(i,subscript))));
            rate(logical(~OffMask(i,subscript)))=rate0(subscript(logical(~OffMask(i,subscript))));
            exchange = rand(1,size(subscript,2)) < rate;
            OffMask(i,subscript(exchange))=~OffMask(i,subscript(exchange));
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