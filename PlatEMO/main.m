function main(varargin)
%main - The interface of PlatEMO.
%
%   main() displays the GUI of PlatEMO.
%
%   main('-Name',Value,'-Name',Value,...) runs one algorithm on a problem
%   with the specified parameter setting.
%
% All the acceptable properties:
%   '-N'            <positive integer>  population size
%   '-M'            <positive integer>  number of objectives
%   '-D'            <positive integer>  number of variables
%	'-algorithm'    <function handle>   algorithm function
%	'-problem'      <function handle>   problem function
%	'-evaluation'   <positive integer>  maximum number of evaluations
%   '-run'          <positive integer>  run number
%   '-save'         <integer>           number of saved populations
%   '-outputFcn'	<function handle>   function invoked after each generation
%
%   Example:
%       main()
%
%   displays the GUI of PlatEMO.
%
%       main('-algorithm',@ARMOEA,'-problem',@DTLZ2,'-N',200,'-M',10)
%
%   runs AR-MOEA on 10-objective DTLZ2 with a population size of 200.
%
%       main('-algorithm',{@KnEA,0.4},'-problem',{@WFG4,6})
%
%   runs KnEA on WFG4, and sets the parameters in KnEA and WFG4.
%
%       for i = 1 : 10
%           main('-algorithm',@RVEA,'-problem',@LSMOP1,'-run',i,'-save',5)
%       end
%
%   runs RVEA on LSMOP1 for 10 times, and each time saves 5 populations to
%   a file in PlatEMO/Data/RVEA.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    cd(fileparts(mfilename('fullpath')));
    addpath(genpath(cd));
    if isempty(varargin)
        if verLessThan('matlab','8.4')
            errordlg('Fail to establish the GUI of PlatEMO, since the version of MATLAB is lower than 8.4 (R2014b). You can run PlatEMO without GUI by invoking main() with parameters.','Error','modal');
        else
            GUI();
        end
    else
        if verLessThan('matlab','7.14')
            error('Fail to execute PlatEMO, since the version of MATLAB is lower than 7.14 (R2012a). Please update the version of your MATLAB software.');
        else
            Global = GLOBAL(varargin{:});
            Global.Start();
        end
    end
end