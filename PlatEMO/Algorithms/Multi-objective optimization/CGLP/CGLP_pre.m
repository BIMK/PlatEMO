function [POP,pop_LCM,pop_DCM,tipe]=CGLP_pre(Problem,hisPop,T,hisPareto,popSize,tipe)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    NP    = size(hisPop{T-1}.X,2);
    Lower = repmat(Problem.lower',1,NP);
    Upper = repmat(Problem.upper',1,NP);
    ran   = Lower'+rand(size(Upper')).*(Upper-Lower)';
    ttyh  = pdist2(hisPop{T-1}.F',hisPop{T-2}.F');
    for ii = 1 : popSize
        [~,mkl(ii)] = min(ttyh(ii,:));
    end
    hisPop{T-2}.X = hisPop{T-2}.X(:,mkl);
    hisPop{T-2}.F = hisPop{T-2}.F(:,mkl);

    %% Correlation analysis
    C41  = hisPop{T-1}.X;  % centroid of time K-1
    C42  = hisPop{T-2}.X;  % centroid of time K-2
    D41  = C41-C42;  % their difference
    C11  = mean(hisPop{T-1}.X',1);  % centroid of time K-1
    C12  = mean(hisPop{T-2}.X',1);  % centroid of time K-2
    D11  = C11'-C12';  % their difference
    D41  = [D41;1:popSize];
    Dd41 = D41;
    num  = cell(3,1);
    te   = [];
    r    = [];
    x    = [Dd41(1:size(C41,1),:),D11]';
    for i = 1 : size(Dd41,2)+1
        x(i,:) = x(i,:)/(x(i,1)+0.0001);
    end
    data = x;
    n    = size(data,1);
    ck   = data(n,:);
    m1   = size(ck,1);
    bj   = data(1:n-1,:);
    m2   = size(bj,1);
    for i = 1 : m1
        for j = 1 : m2
            te(j,:) = bj(j,:)-ck(i,:);
        end
        jc1 = min(min(abs(te')));
        jc2 = max(max(abs(te')));
        rho = 0.5;
        ksi = (jc1+rho*jc2)./(abs(te)+rho*jc2);
        rt  = sum(ksi')/size(ksi,2);
        r(i,:) = rt;
    end
    [~,rind] = sort(r,'descend');
    
    %% Divide three groups
    num{1} = Dd41(end,rind(1:6/10*popSize));
    num{2} = Dd41(end,rind(6/10*popSize+1:9/10*popSize));
    num{3} = Dd41(end,rind(9/10*popSize+1:end));
    pre_solution = zeros(popSize,size(Lower,1));
    
    %% Operator for high correlation group
    opp  = 1;
    Cw11 = mean(hisPop{T-1}.X(:,num{opp})',1);  % centroid of time K-1
    Cw12 = mean(hisPop{T-2}.X(:,num{opp})',1);  % centroid of time K-2
    Dw11 = Cw11'-Cw12';  % their difference
    are  = repmat(Dw11',size(num{opp},2),1);
    pre_solution(num{opp},:) = hisPop{T-1}.X(:,num{opp})'+are;
    
    %% Operator for mid correlation group
    pop_LCM = [];
    pop_DCM = [];
    [pre_solution(num{2},:),pop_LCM,pop_DCM,tipe] = DLCM(hisPop{T-1}.X(:,num{2}), hisPop{T-2}.X(:,num{2}),pop_LCM,pop_DCM,tipe);
    
    %% Operator for low correlation group
    if size(hisPareto{T-1},1) > size(num{end},2)
        pre_solution(num{end},:) = hisPareto{T-1}(1:size(num{end},2),:);
    else
        pre_solution(num{end}(1:size(hisPareto{T-1},1)),:) = hisPareto{T-1};
    end
    POP = pre_solution(1:popSize,:);
    POP(POP<Lower'|POP>Upper') = ran(POP<Lower'|POP>Upper');
end