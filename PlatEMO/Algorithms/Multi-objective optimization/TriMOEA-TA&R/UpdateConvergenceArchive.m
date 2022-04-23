function [AC,Rank,fS] = UpdateConvergenceArchive(AC,Population,NC,Z,Xic,sigma_niche,Problem)
% Update the Convergence Archive

%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    Population = [Population,AC];
    PopObj = Population.objs;
    PopDec = Population.decs;
    N = size(PopObj,1);   
    Rank   = inf(1,N);      % Rank of each solution
    nrank  = 1;             % Current rank
    
    %% Normalize
    PopObj = PopObj - repmat(Z,N,1);
    PopDec = (PopDec - repmat(Problem.lower,N,1))./repmat(Problem.upper-Problem.lower,N,1);
    
    %% Convergence indicator
    fS = mean(PopObj');
    
    %% Calculate distance between every two solutions in the IC decsion subspace
    d = pdist2(PopDec(:,Xic),PopDec(:,Xic),'chebychev');
       
    %% Rank
    Choose = false(1,N);
    Q = true(1,N);
    Q1 = false(1,N);
    while sum(Choose) < NC
        if sum(Q) == 0
           Q = Q1;
           Q1 = false(1,N);
           nrank = nrank+1;
        end
        % Choose x with min FS 
        temp1 = fS == min(fS(Q));
        xmin = find(and(temp1,Q));
        xmin = xmin(1);
        Rank(xmin) = nrank;
        Choose(xmin) = true;
        Q(xmin) = false;        
        % Delete solution near x_min
        temp3=d(xmin,:);
        temp2=temp3<sigma_niche;
        Delete = and(temp2,Q);
        Q(Delete) = false;
        Q1(Delete) = true;      
    end
    AC = Population(Choose);
    Rank = Rank(Choose);
    fS = fS(Choose);
end