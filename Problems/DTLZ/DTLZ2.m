function varargout = DTLZ2(Operation,Global,input)
% <problem> <DTLZ>
% Scalable Test Problems for Evolutionary Multi-Objective Optimization
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        %% Initialize the parameter setting, and randomly generate a population
        case 'init'
            % Set the default number of objectives
            Global.M        = 3;
            % Set the default number of decision variables
            Global.D        = Global.M + 9;
            % Set the lower bound of each decision variable
            Global.lower    = zeros(1,Global.D);
            % Set the upper bound of each decision variable
            Global.upper    = ones(1,Global.D);
            % Set the default operator for this problem
            Global.operator = @EAreal;
            % Randomly generate a number of solutions' decision variables
            PopDec    = rand(input,Global.D);
            % Return the set of decision variables
            varargout = {PopDec};
        %% Calculate the objective values of a population
        case 'value'
            % The set of decision variables
            PopDec = input;
            % The number of objectives
            M      = Global.M;
            % Calculate the objective values according to the decision
            % variables, note that here the objective values of multiple
            % solutions are calculated at the same time
            g      = sum((PopDec(:,M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            % Calculate the constraint values
            PopCon = [];
            % Return the decision variables, objective values, and
            % constraint values
            varargout = {input,PopObj,PopCon};
        %% Generate reference points on the true Pareto front
        case 'PF'
            % Generate a set of reference points on the true Pareto front
            f = UniformPoint(input,Global.M);
            f = f./repmat(sqrt(sum(f.^2,2)),1,Global.M);
            % Return the reference points
            varargout = {f};
    end
end