function [Population1,Mask1,Dec1,FrontNo1,SV1,Population2,Mask2,Dec2,FrontNo2,SV2] = Divide(Population1,Mask1,Dec1,FrontNo1,row,col)
% The divide operation of HHC-MMEA
%The Divisionflag() function returns the sequence number to perform this split

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % subpopulation2
    Population2 = Population1(col);
    Dec2 = Dec1(col,:);
    Mask2 = Mask1(col,:);
    FrontNo2 = FrontNo1(col);
    mask = Mask2(FrontNo2==1,:);
    SV2 = sum(mask,1)/size(mask,1);

    %subpopulation1
    Population1 = Population1(row);
    Dec1 = Dec1(row,:);
    Mask1= Mask1(row,:);
    FrontNo1 = FrontNo1(row);
    mask = Mask1(FrontNo1==1,:);
    SV1 = sum(mask,1)/size(mask,1);
end