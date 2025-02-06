function NewDec = Match(Dec,Mask,Problem)
% Matching Dec and Mask

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Nomarlize
    NDec   = (Dec - Problem.lower) ./ (Problem.upper - Problem.lower);
    NewDec = zeros(size(Mask,1),Problem.D);
    index  = 1 : 1 : size(Mask,1);
    
    % Compute the similarity between the normalized Dec and Mask and match
    while ~isempty(index)
        k      = randi(length(index));
        Cosine = 1 - pdist2(Mask(index(k),:),NDec,'cosine');
        [~,H]  = max(Cosine,[],2);
        NewDec(index(k),:) = NDec(H,:);
        index(k)  = [];
        NDec(H,:) = [];
    end
    
    % Repair
    NewDec = (NewDec .* (Problem.upper - Problem.lower)) + Problem.lower;
end