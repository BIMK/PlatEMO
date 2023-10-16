function [LDis,B] = Initialize_SOM(S,D,H)
% Initialize the SOM

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Chao He

   % Position of each neuron
    D = arrayfun(@(S)1:S,D,'UniformOutput',false);
    eval(sprintf('[%s]=ndgrid(D{:});',sprintf('c%d,',1:length(D))))
    eval(sprintf('Z=[%s];',sprintf('c%d(:),',1:length(D))))
    % Distance between each two neurons in latent space
    LDis = pdist2(Z,Z);
    % H nearest neurons of each neuron in latent space
    [~,B] = sort(LDis,2);
    B     = B(:,2:min(H+1,end));
end