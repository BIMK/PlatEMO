function index = TournamentSelection(K,N,varargin)
%TournamentSelection - Tournament selection.
%
%   P = TournamentSelection(K,N,fitness1,fitness2,...) returns the indices
%   of N solutions by K-tournament selection based on their fitness values.
%   In each selection, the candidate having the minimum fitness1 value will
%   be selected; if multiple candidates have the same minimum value of
%   fitness1, then the one with the smallest fitness2 value is selected,
%   and so on.
%
%   Example:
%       P = TournamentSelection(2,100,FrontNo)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    varargin    = cellfun(@(S)reshape(S,[],1),varargin,'UniformOutput',false);
    [Fit,~,Loc] = unique([varargin{:}],'rows');
    [~,rank]    = sortrows(Fit);
    [~,rank]    = sort(rank);
    Parents     = randi(length(varargin{1}),K,N);
    [~,best]    = min(rank(Loc(Parents)),[],1);
    index       = Parents(best+(0:N-1)*K);
end