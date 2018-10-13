function CSEA(Global)
% <algorithm> <A-G>
% A Classification Based Surrogate-Assisted Evolutionary Algorithm for
% Expensive Many-Objective Optimization
% k    ---    6 --- Number of reference solutions
% gmax --- 3000 --- Number of solutions evaluated by surrogate model

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
    [k,gmax] = Global.ParameterSet(6,3000);
    
    %% Initalize the population by Latin hypercube sampling
    N          = 11*Global.D-1;
    PopDec     = lhsamp(N,Global.D);
    Population = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
    Arc        = Population;
    
    %% Optimization
	while Global.NotTermination(Arc)
        % Select reference solutions and preprocess the data
        Ref    = RefSelect(Population,k);
        Input  = Population.decs;  
        Output = GetOutput(Population.objs,Ref.objs); 
        rr     = sum(Output)/length(Output);
        tr     = min(rr,1-rr)*0.5;
        [TrainIn,TrainOut,TestIn,TestOut] = DataProcess(Input,Output);
        
        % Construct and update the FNN
        net = NN(ceil(Global.D*2),800);
        net.train(TrainIn,TrainOut);
        
        % Error rates calculation
        TestPre   = net.predict(TestIn);
        IndexGood = TestOut==1;
        p0 = sum(abs((TestOut(IndexGood)-TestPre(IndexGood))))/sum(IndexGood);
        p1 = sum(abs((TestOut(~IndexGood)-TestPre(~IndexGood))))/sum(~IndexGood);
        
        % Surrogate-assisted selection and update the population
        Next = SurrogateAssistedSelection(Global,net,p0,p1,Ref,Population.decs,gmax,tr);
        if ~isempty(Next)
            Arc = [Arc,INDIVIDUAL(Next)];
        end
        Population = RefSelect(Arc,Global.N);
	end
end