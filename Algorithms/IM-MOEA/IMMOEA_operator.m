function Offspring = IMMOEA_operator(Global,Population)
% <operator> <real>
% The Gaussian process based reproduction
% L --- 3 --- Size of random group

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is modified from the code in
% http://www.soft-computing.de/jin-pub_year.html

    L = Global.ParameterSet(3);
    PopDec = Population.decs;
	PopObj = Population.objs;
    [N,D]  = size(PopDec);
    
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
    randDec = rand(N,D).*(Upper-Lower) + Lower;
    invalid = OffDec<Lower | OffDec>Upper;
    OffDec(invalid) = randDec(invalid);

    %% Polynomial mutation
    [proM,disM] = deal(1,20);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    
    Offspring = INDIVIDUAL(OffDec);
end