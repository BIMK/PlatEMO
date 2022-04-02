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
%   'save'       	<integer>           number of saved populations
%   'outputFcn'     <function handle>   function called after each iteration
%   'encoding'      <string>            encoding scheme of variables
%   'lower'         <vector>            lower bounds of variables
%   'upper'         <vector>            upper bounds of variables
%   'decFcn'        <function handle>   function of variable repair
%   'objFcn'        <function handle>   function of objective calculation
%   'conFcn'        <function handle>   function of constraint calculation
%
%   Example:
%
%       platemo()
%
%   displays the GUI of PlatEMO.
%
%       platemo('algorithm',@GA,'problem',@SOP_F1,'-N',50)
%
%   runs GA with a population size of 50 on SOP_F1.
%
%       platemo('algorithm',{@KnEA,0.4},'problem',{@WFG4,6})
%
%   runs KnEA on WFG4 and sets the parameters in KnEA and WFG4.
%
%       for i = 1 : 10
%           platemo('algorithm',@MOEAD,'problem',@ZDT1,'save',5)
%       end
%
%   runs MOEA/D on ZDT1 for 10 times, where 5 populations are saved to a
%   file in PlatEMO/Data/MOEAD in each time.
%
%       platemo('algorithm',@CCMO,'objFcn',@CalObj,'conFcn',@CalCon)
%
%   runs CCMO on a problem whose objective values are calculated by
%   CalObj() and constraint violations are calculated by CalCon().

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
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
                Algorithm = ALG(input{:},'save',-1);
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

function [name,input] = getSetting(setting,Pro)
    isStr = find(cellfun(@ischar,setting(1:end-1))&~cellfun(@isempty,setting(2:end)));
    input = {};
    if nargin > 1
        index = isStr(find(strcmp(setting(isStr),'algorithm'),1)) + 1;
        if isempty(index)
            names = {@BSPGA,@GA,@SACOSO,@GA;@PMMOEA,@NSGAIII,@KRVEA,@NSGAIII;@RVEA,@RVEA,@CSEA,@RVEA};
            name  = names{find([Pro.M<2,Pro.M<4,1],1),find([ismember({'binary','permutation'},Pro.encoding),Pro.maxFE<=1000&Pro.D<=10,1],1)};
        elseif iscell(setting{index})
            name  = setting{index}{1};
            input = {'parameter',setting{index}(2:end)};
        else
            name = setting{index};
        end
        keys = {'save','outputFcn'};
    else
        index = isStr(find(strcmp(setting(isStr),'problem'),1)) + 1;
        if isempty(index)
            name = @UserProblem;
            keys = {'N','D','maxFE','encoding','lower','upper','initFcn','decFcn','objFcn','conFcn','parameter'};
        elseif iscell(setting{index})
            name  = setting{index}{1};
            input = {'parameter',setting{index}(2:end)};
            keys  = {'N','M','D','maxFE'};
        else
            name = setting{index};
            keys = {'N','M','D','maxFE'};
        end
    end
    index = isStr(ismember(setting(isStr),keys));
    input = [input,setting(reshape([index;index+1],1,[]))];
end