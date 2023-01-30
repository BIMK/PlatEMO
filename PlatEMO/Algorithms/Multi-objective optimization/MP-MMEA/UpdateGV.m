function gv = UpdateGV(gv,Mask,FrontNo)
% Update the guiding vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Mask = Mask(FrontNo==1,:);
    k    = ceil(0.1*size(Mask,1));
    v    = zeros(1,size(Mask,2));
    for i = 1 : size(Mask,1)
        rand_Mask = Mask(i,:);          
        dis       = pdist2(rand_Mask,Mask,'hamming');
        [~,index] = sort(dis);
        knearest_Mask = Mask(index(1:k),:);
        rand_Mask(rand_Mask==1) = 0.5 ;
        knearest_Mask(knearest_Mask==1) = 0.5 ;
        v = v + sum(repmat(rand_Mask,k,1)+knearest_Mask,1)/(size(Mask,1)*k);
    end
    if all(gv==0)
        gv = v;
    else
        gv = 0.9*gv + 0.1*v;
    end
end