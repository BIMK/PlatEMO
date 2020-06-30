function varargout = GridDensity(div,varargin)
% Calculate the number of solutions in the grid of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the grid location of each solution
    PopObj = cat(1,varargin{:});
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    d      = (fmax-fmin)/div;
    GLoc   = floor((PopObj-repmat(fmin,size(PopObj,1),1))./repmat(d,size(PopObj,1),1));
    GLoc(GLoc>=div)   = div - 1;
    GLoc(isnan(GLoc)) = 0;
    
    %% Calculate the number of solutions in the grid of each solution
    [~,~,Site] = unique(GLoc,'rows');
    CrowdG     = hist(Site,1:max(Site));
    Crowd      = CrowdG(Site);
    varargout  = mat2cell(Crowd,1,cellfun(@(S)size(S,1),varargin));
end