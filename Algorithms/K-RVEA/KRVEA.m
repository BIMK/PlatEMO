function KRVEA(Global)
% <algorithm> <H-N>
% A Surrogate-assisted Reference Vector Guided Evolutionary Algorithm for
% Computationally Expensive Many-objective Optimization
% alpha ---  2 --- The parameter controlling the rate of change of penalty
% wmax  --- 20 --- Number of generations before updating Kriging models
% mu    ---  5 --- Number of re-evaluated solutions at each generation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Parameter setting
    [alpha,wmax,mu] = Global.ParameterSet(2,20,5);

    %% Generate the reference points and population
    [V0,Global.N] = UniformPoint(Global.N,Global.M);
	V     = V0;
    NI    = 11*Global.D-1;
    P     = lhsamp(NI,Global.D);
    A2    = INDIVIDUAL(repmat(Global.upper-Global.lower,NI,1).*P+repmat(Global.lower,NI,1));
    A1    = A2;  
    THETA = 5.*ones(Global.M,Global.D);
    Model = cell(1,Global.M);
    
    %% Optimization
    while Global.NotTermination(A2)
        % Refresh the model and generate promising solutions
        A1Dec = A1.decs;
        A1Obj = A1.objs;
        for i = 1 : Global.M
            % The parameter 'regpoly1' refers to one-order polynomial
            % function, and 'regpoly0' refers to constant function. The
            % former function has better fitting performance but lower
            % efficiency than the latter one
            dmodel     = dacefit(A1Dec,A1Obj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Global.D),100.*ones(1,Global.D));
            Model{i}   = dmodel;
            THETA(i,:) = dmodel.theta;
        end
        PopDec = A1Dec;
        w      = 1;
        while w <= wmax
            drawnow();
            OffspringDec = Global.VariationDec(PopDec);
            PopDec = [PopDec;OffspringDec];
            [N,~]  = size(PopDec);
            PopObj = zeros(N,Global.M);
            MSE    = zeros(N,Global.M);
            for i = 1: N
                for j = 1 : Global.M
                    [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Model{j});
                end
            end
            index  = KEnvironmentalSelection(PopObj,V,(w/wmax)^alpha);
            PopDec = PopDec(index,:);
            PopObj = PopObj(index,:);
            % Adapt referece vectors
            if ~mod(w,ceil(wmax*0.1))
                V(1:Global.N,:) = V0.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V0,1),1);
            end
            w = w + 1; 
        end
        
        % Select mu solutions for re-evaluation
        [NumVf,~] = NoActive(A1Obj,V0);
        PopNew    = KrigingSelect(PopDec,PopObj,MSE,V,V0,NumVf,0.05*Global.N,mu,(w/wmax)^alpha); 
        New       = INDIVIDUAL(PopNew);
        A2        = [A2,New];
        A1        = UpdataArchive(A1,New,V,mu,NI); 
    end
end