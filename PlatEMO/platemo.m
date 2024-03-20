function varargout = platemo(varargin)
%platemo - The main function of PlatEMO.
%
%   platemo() displays the GUI of PlatEMO.
%
%   platemo('Name',Value,'Name',Value,...) runs an algorithm on a problem
%   with specified parameter settings.
%
% All the acceptable names and values are:
%	'algorithm'     <function handle>	an algorithm
%	'problem'       <function handle>	a problem
%   'N'             <positive integer>  population size
%   'M'             <positive integer>  number of objectives
%   'D'             <positive integer>  number of variables
%	'maxFE'         <positive integer>  maximum number of function evaluations
%   'maxRuntime'    <positive integer>  maximum runtime (in second)
%   'save'       	<integer>           number of saved populations
%   'run'           <integer>           current run number
%   'metName'       <string>            names of metrics to calculate
%   'outputFcn'     <function handle>   function called after each iteration
%   'encoding'      <string>            encoding scheme of each decision variable (1.real 2.integer 3.label 4.binary 5.permutation)
%   'lower'         <vector>            lower bound of each decision variable
%   'upper'         <vector>            upper bound of each decision variable
%   'initFcn'       <function handle>   function for initializing solutions
%   'evalFcn'       <function handle>   function for evaluating solutions
%   'decFcn'        <function handle>   function for repairing invalid solutions
%   'objFcn'        <function handle>   objective functions
%   'conFcn'        <function handle>   constraint functions
%   'objGradFcn'    <function handle>   function for calculating the gradients of objectives
%   'conGradFcn'    <function handle>   function for calculating the gradients of constraints
%
%   Example:
%
%       platemo()
%
%   displays the GUI of PlatEMO.
%
%       platemo('algorithm',@GA,'problem',@SOP_F1,'N',50,'maxFE',20000)
%
%   runs GA with a population size of 50 on SOP_F1 for 20000 evaluations.
%
%       platemo('algorithm',@PSO,'problem',@SOP_F1,'N',100,'maxRuntime',3)
%
%   runs PSO with a population size of 100 on SOP_F1 for 3 seconds.
%
%       platemo('algorithm',{@KnEA,0.4},'problem',{@WFG4,6},'M',5)
%
%   runs KnEA on 5-objective WFG4 and sets the parameters in KnEA and WFG4.
%
%       for i = 1 : 10
%           platemo('algorithm',@MOEAD,'problem',@ZDT1,'save',5,'run',i)
%       end
%
%   runs MOEA/D on ZDT1 for 10 times, where 5 populations are saved to a
%   distinct file in PlatEMO/Data/MOEAD each time.
%
%       platemo('algorithm',@CCMO,'objFcn',{@Obj1,@Obj2},'conFcn',@Con)
%
%   runs CCMO on a problem whose two objective values are calculated by
%   Obj1() and Obj2() and one constraint violation is calculated by Con().
%
%       platemo('algorithm',@SparseEA,'evalFcn',@Evaluation)
%
%   runs SparseEA on a problem whose objective values and constraint
%   violations are all calculated by a single function Evaluation().

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    cd(fileparts(mfilename('fullpath')));
    addpath(genpath(cd));
    if isempty(varargin)
        if verLessThan('matlab','9.9')
            errordlg('Fail to create the GUI of PlatEMO since the version for MATLAB is lower than R2020b. You can use PlatEMO without GUI by calling platemo() with parameters.','Error','modal');
        else
            try
                GUI();
            catch err
                errordlg('Fail to create the GUI, please make sure all the folders of PlatEMO have been added to search path.','Error','modal');
                rethrow(err);
            end
        end
    else
        if verLessThan('matlab','9.4')
            error('Fail to use PlatEMO since the version for MATLAB is lower than R2018a. Please update your MATLAB software.');
        else
            [PRO,input] = getSetting(varargin);
            Problem     = PRO(input{:});
            [ALG,input] = getSetting(varargin,Problem);
            if nargout > 0
                Algorithm = ALG(input{:},'save',0);
            else
                Algorithm = ALG(input{:});
            end
            Algorithm.Solve(Problem);
            if nargout > 0
                P = Algorithm.result{end};
                varargout = {P.decs,P.objs,P.cons};
            end
        end
    end
end

function [name,Setting] = getSetting(Setting,Pro)
    isStr = find(cellfun(@ischar,Setting(1:end-1))&~cellfun(@isempty,Setting(2:end)));
    if nargin > 1
        index = isStr(find(strcmp(Setting(isStr),'algorithm'),1)) + 1;
        if isempty(index)
            names = {@BSPGA,@GA,@SACOSO,@GA;@PMMOEA,@NSGAIII,@KRVEA,@NSGAIII;@RVEA,@RVEA,@CSEA,@RVEA};
            name  = names{find([Pro.M<2,Pro.M<4,1],1),find([all(Pro.encoding==4),any(Pro.encoding>2),Pro.maxFE<=1000&Pro.D<=10,1],1)};
        elseif iscell(Setting{index})
            name    = Setting{index}{1};
            Setting = [Setting,{'parameter'},{Setting{index}(2:end)}];
        else
            name = Setting{index};
        end
    else
        index = isStr(find(strcmp(Setting(isStr),'problem'),1)) + 1;
        if isempty(index)
            name = @UserProblem;
        elseif iscell(Setting{index})
            name    = Setting{index}{1};
            Setting = [Setting,{'parameter'},{Setting{index}(2:end)}];
        else
            name = Setting{index};
        end
    end
end