function dim_base = DimBase(Mask_dim,dim,num)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    non_zero = zeros(1,size(Mask_dim,1));  
    for i = 1 : size(Mask_dim,1)
        non_zero(i) = size(find(Mask_dim(i,:)),2);
    end
    index = find(non_zero==num);
    n     = size(index,2); 
    dim_b = zeros(n,num); 
    for i = 1 : n
        dim_b(i,:) = dim(logical(Mask_dim(index(i),:)));
    end 
    dim_base = unique(dim_b)'; 
end