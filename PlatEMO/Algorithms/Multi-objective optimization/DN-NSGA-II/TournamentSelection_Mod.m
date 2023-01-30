function index = TournamentSelection_Mod(K,N,PopDec,varargin)
% Tournament selection of DN-NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    varargin = cellfun(@(S)reshape(S,length(varargin{1}),1),varargin,'UniformOutput',false);
    [~,rank] = sortrows([varargin{:}]);
    [~,rank] = sort(rank);    
    Parents  = randi(length(varargin{1}),K,N);
    for i = 1 : N
        [~,min_index] = min(pdist2(PopDec(Parents(1,i),:),PopDec(Parents(2:K,i),:)));
        Parents(2,i)  = Parents(min_index+1,i);
    end 
    Parents  = Parents(1:2,:);
    [~,best] = min(rank(Parents),[],1);
    index    = Parents(best+(0:N-1)*2);
end