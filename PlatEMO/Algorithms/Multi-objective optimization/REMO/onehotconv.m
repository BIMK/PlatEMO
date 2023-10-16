function varargout = onehotconv(varargin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if varargin{2}== 1
        %% conv onehot
        l        = varargin{1};      
        l_onehot = zeros(size(l,1),3);
     
        l_onehot(l == 1 ,1) = 1;
        l_onehot(l == 0,2)  = 1;
        l_onehot(l == -1,3) = 1;
        
        varargout = {l_onehot};
        
    elseif varargin{2} == 2
        %% deconv onehot
        onehot_l = varargin{1};
        res_l    = zeros(size(onehot_l,1),1);

        [~,maxind] = max(onehot_l,[],2);

        res_l(maxind==1) = 1;
        res_l(maxind==3) = -1;

        varargout = {res_l};
    end
end