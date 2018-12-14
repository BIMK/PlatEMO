function Archive = UpdateArchive(Archive,Population,popsize)
% Update the archive in MSOPS-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Combine the archive with the population
    Archive = [Archive,Population];
    ArchObj = Archive.objs;
    N       = length(Archive);

    %% Update the archive by weighted min-max metric
    % Calculate the weighted min-max metric between each two solutions
    WMM = CalMetric(ArchObj,ArchObj);
    WMM_diagonal = WMM(logical(eye(N)))';
    % Delete the solutions which do not have the lowest metric value than
    % others according to its own weight vector
    Remain = true(1,N);
    for i = 1 : N
        if Remain(i)
            if WMM(i,i) > min(WMM(Remain,i))
                Remain(i) = false;
            else
                Remain(WMM(i,:)<WMM_diagonal) = false;
            end
        end
    end
    Archive = Archive(Remain);
    % If the archive contains too many solutions, randomly delete some. The
    % original algorithm does not limit the size of archive, so that the
    % size of archive will increase unrestrainedly
    if length(Archive) > 10*popsize
        Archive = Archive(randperm(length(Archive),5*popsize));
    end
end