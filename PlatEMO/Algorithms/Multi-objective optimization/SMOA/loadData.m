function Population = loadData(Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    dir = fullfile(pwd,'Algorithms','Multi-objective optimization','SMOA');
    [file,path] = uigetfile('*.*','Select a Supervised Data File',dir);

    if endsWith(file,'.mat') % PlarEMO's default save file
        load(fullfile(path, file),'result');
        Population = result{end,2};
    else % Decision variables save file
        Dec = load(fullfile(path, file));
        Population = Problem.Evaluation(Dec);
    end
end