function Offspring = Operator(Population,L)
% The Gaussian process based reproduction

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is modified from the code in
% http://www.soft-computing.de/jin-pub_year.html

    %% Parameter setting
    if nargin < 2
        L = 3;
    end
    PopDec = Population.decs;
	PopObj = Population.objs;
    [N,D]  = size(PopDec);
    Global = GLOBAL.GetObj();
    
    %% Gaussian process based reproduction
    if length(Population) < 2*Global.M
        OffDec = PopDec;
    else
        OffDec = [];
        fmin   = 1.5*min(PopObj,[],1) - 0.5*max(PopObj,[],1);
        fmax   = 1.5*max(PopObj,[],1) - 0.5*min(PopObj,[],1);
        % Train one groups of GP models for each objective
        for m = 1 : Global.M
            parents = randperm(N,floor(N/Global.M));
            offDec  = PopDec(parents,:);
            for d = randperm(D,L)
                % Gaussian Process
                try
                    [ymu,ys2] = gp(struct('mean',[],'cov',[],'lik',log(0.01)),...
                                   @infExact,@meanZero,@covLIN,@likGauss,...
                                   PopObj(parents,m),PopDec(parents,d),...
                                   linspace(fmin(m),fmax(m),size(offDec,1))');
                    offDec(:,d) = ymu + rand*sqrt(ys2).*randn(size(ys2));
                catch
                end
            end
            OffDec = [OffDec;offDec];
        end
    end
    
    %% Convert invalid values to random values
    [N,D]   = size(OffDec);
    Lower   = repmat(Global.lower,N,1);
    Upper   = repmat(Global.upper,N,1);
    randDec = unifrnd(Lower,Upper);
    invalid = OffDec<Lower | OffDec>Upper;
    OffDec(invalid) = randDec(invalid);

    %% Polynomial mutation
    [proM,disM] = deal(1,20);
    Site   = rand(N,D) < proM/D;
    mu     = rand(N,D);
    temp   = Site & mu<=0.5;
    OffDec = min(max(OffDec,Lower),Upper);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = INDIVIDUAL(OffDec);
end