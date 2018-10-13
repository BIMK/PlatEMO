function main(varargin)
%main - The interface of PlatEMO
%
%   main() runs the GUI of PlatEMO.
%
%   main('-Name',Value,'-Name',Value,...) runs one algorithm on a problem
%   with the specified parameter setting.
%
% All the acceptable properties:
%   '-N'            <positive integer>  population size
%   '-M'            <positive integer>  number of objectives
%   '-D'            <positive integer>  number of variables
%	'-evaluation'   <positive integer>  maximum number of evaluations
%	'-algorithm'    <function handle>   algorithm function
%	'-problem'      <function handle>   problem function
%	'-operator'     <function handle>   operator function
%   '-mode'         <positive integer>  run mode (1.show result 2.save result 3.run outputFcn)
%   '-run'          <positive integer>  run No.
%   '-outputFcn'	<function handle>   function invoked after each generation when mode = 3
%	'-X_parameter'  <cell>              the parameter values of function X
%
%   Example:
%       main()
%       main('-algorithm',@NSGAII,'-problem',@DTLZ2,'-N',100,'-M',2)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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