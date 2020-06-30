function A = UpdateArchive(A,S,K)
% Update the external archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Combine A and S and normalize the objective values
    Combine = [A,S];
    fmin    = repmat(min(A.objs,[],1),length(Combine),1);
    fmax    = repmat(max(A.objs,[],1),length(Combine),1);
    PopObj  = (Combine.objs-fmin)./(fmax-fmin);
    [N,M]   = size(PopObj);
    
    %% Calculate the shifted distance between each two solutions
    sde = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            sde(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    
    %% Calculate Cv
    dis = sqrt(sum(PopObj.^2,2));
    % Use 1-dis instead of 1-dis/sqrt(M)
    Cv  = 1 - dis;
    
    %% Calculate d1 and d2
    Cosine = 1 - pdist2(PopObj,ones(1,M),'cosine');
    d1     = dis.*Cosine;
    d2     = dis.*sqrt(1-Cosine.^2);

    %% Insert each solution in S into the archive A
    Choose = 1 : length(A); % Indices of selected solutions in [A,S]
    for i = 1 : length(S)
        % Check the dominance relations between S(i) and solutions in A
        mark = false(1,length(Choose)+1);
        for j = 1 : length(Choose)
            flag = any(S(i).obj<Combine(Choose(j)).obj,2) - any(S(i).obj>Combine(Choose(j)).obj,2);
            if flag == 1
                mark(j) = true;
            elseif flag == -1
                mark(end) = true;
                break;
            end
        end
        Choose(mark(1:end-1)) = [];
        if ~mark(end)
            % Insert S(i) into A
            Choose = [Choose,length(A)+i];
            if length(Choose) > K
                % Delete the one in A with the lowest BFE value
                [~,worst]     = min(CalBFE(sde(Choose,Choose),Cv(Choose),d1(Choose),d2(Choose)));
                Choose(worst) = [];
            end
        end
    end
    
    %% Sort the solutions in A according to their BFE values
    A        = Combine(Choose);
    [~,rank] = sort(CalBFE(sde(Choose,Choose),Cv(Choose),d1(Choose),d2(Choose)),'descend');
    A        = A(rank);
end

function BFE = CalBFE(sde,Cv,d1,d2)
% Calculate the BFE value of each solution in the population

% This function is modified from the code in
% http://security.szu.edu.cn/people.aspx?p=QiuzhenLin
% Where the parameters are a little different from the ones in the paper

    %% Calculate Cd
    SDE = min(sde,[],2);
    Cd  = (SDE-min(SDE))./(max(SDE)-min(SDE));
    
    %% Determine the value of alpha and beta of each solution
    alpha   = zeros(length(Cv),1);
    beta    = zeros(length(Cv),1);
    meanCd  = mean(Cd);
    meanCv  = mean(Cv);
    meand1  = mean(d1);
    meand2  = mean(d2);
    case111 = Cv >  meanCv & d1 <= meand1 & Cd <= meanCd;
    case112 = Cv >  meanCv & d1 <= meand1 & Cd >  meanCd;
    case121 = Cv >  meanCv & d1 >  meand1 & Cd <= meanCd;
    case122 = Cv >  meanCv & d1 >  meand1 & Cd >  meanCd;
    case211 = Cv <= meanCv & d1 <= meand1 & d2 >  meand2 & Cd <= meanCd;
    case212 = Cv <= meanCv & d1 <= meand1 & d2 >  meand2 & Cd >  meanCd;
    case221 = Cv <= meanCv &(d1 >  meand1 | d2 <= meand2)& Cd <= meanCd;
    case222 = Cv <= meanCv &(d1 >  meand1 | d2 <= meand2)& Cd >  meanCd;
    alpha(case111) = rand(sum(case111),1)*0.3+0.8; beta(case111) = 1;
    alpha(case112) = 1;   beta(case112) = 1;
    alpha(case121) = 0.6; beta(case121) = 1;
    alpha(case122) = 0.9; beta(case122) = 1;
    alpha(case211) = rand(sum(case211),1)*0.3+0.8; beta(case211) = rand(sum(case211),1)*0.3+0.8;
    alpha(case212) = 1;   beta(case212) = 1;
    alpha(case221) = 0.2; beta(case221) = 0.2;
    alpha(case222) = 1;   beta(case222) = 0.2;

    %% The BFE value of each solution
    BFE = alpha.*Cd + beta.*Cv;
end