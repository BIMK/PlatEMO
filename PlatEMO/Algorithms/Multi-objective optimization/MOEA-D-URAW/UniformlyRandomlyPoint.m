function [W1,N] = UniformlyRandomlyPoint(N,M)
%UniformlyRandomlyPoint - Generate a set of uniform randomly distributed points on
%the unit hyperplane
%
%   [W,N] = UniformlyRandomlyPoint(N,M) returns N uniform randomly distributed
%   points with M objectives.
%
%   Example:
%       [W,N] = UniformlyRandomlyPoint(275,10)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias
	
    W1=eye(M,M);
	W1=[W1;ones(1,M)/M];
    
	W2=rand(5000,M);
	W2 = W2./repmat(sum(W2,2),1,size(W2,2));
	
	while size(W1,1) < N
		index = find_index_with_largest_distance (W1,W2);
		W1(size(W1,1)+1,:)=W2(index,:);
		W2(index,:)=[];
	end	

    W1 = max(W1,1e-6);
    N = size(W1,1);

end

function index = find_index_with_largest_distance (W1,W2)
    Distance = pdist2(W2,W1);
    Temp     = sort(Distance,2);
    [~,Rank] = sortrows(Temp);
    index=Rank(length(Rank));
end