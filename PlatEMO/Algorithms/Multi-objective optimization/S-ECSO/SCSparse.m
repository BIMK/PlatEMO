function Population = SCSparse(Problem,Population,L)
% The sparse operator of S-ECSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    x = Population.decs;
    
    %% Make the population sparse
    [row_x,col_x] = size(x);
    new_x = [];
    for i = 1 : row_x
        new_x(i,:) = zeros(1,col_x);
        flag_T = L < abs(x(i,:));
        new_x(i,flag_T) = x(i,flag_T)- L(flag_T).*sign(x(i,flag_T));
    end

    %% Restrict the range
    x_maxmin = (Problem.upper-Problem.lower)';
    for irange = 1:row_x %for every solution
        Upper_flag = Problem.upper<new_x(irange,:);
        Upper_flag_T = sum(Upper_flag);
        if Upper_flag_T > 0
            Upper_index = find(Upper_flag == 1);
            new_x(irange,Upper_index) = rand(1,size(Upper_index,2)).*x_maxmin(Upper_index);
        end

        Low_flag = Problem.lower > new_x(irange,:);
        Low_flag_T = sum(Low_flag);
        if Low_flag_T > 0
            Low_index = find(Low_flag == 1);
            new_x(irange,Low_index) = rand(1,size(Low_index,2)).*x_maxmin(Low_index);
        end
    end

    %% Calculate obj
    new_x = Problem.Evaluation(new_x);
    XX    = [Population, new_x];
    
    %% Save the better solutions
    N = length(new_x);
    [FrontNo,~] = NDSort([XX.objs],Inf);
    for i = 1 : N
        if FrontNo(i+N) <= FrontNo(i)
            Population(i) =  XX(i+N);
        end
    end
end