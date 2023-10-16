function [Archive,ArcW,ArcSP] = UpdateArchive(Archive,ArcObjX,SOI,SOIObjX,alpha,z,W,ArcW,ArcSP)
% Update the archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Update the archive
    Tr      = max(sum(abs(SOIObjX-z),2));
    Remain  = sum(abs(ArcObjX-z),2) <= alpha*Tr;
    Archive = [Archive(Remain),SOI];
    
    %% Update other information
    [~,SOIW]    = min(pdist2(SOIObjX-z,W,'cosine'),[],2);
    [~,ArcWNew] = min(pdist2(ArcObjX-z,W,'cosine'),[],2);
    ArcW        = arrayfun(@(i)[ArcW{i},ArcWNew(i)],(1:length(ArcW))','UniformOutput',false);
    ArcW        = [ArcW(Remain);num2cell(SOIW)];
    ArcSPNew    = sum(ArcObjX,2);
    ArcSP       = arrayfun(@(i)[ArcSP{i},ArcSPNew(i)],(1:length(ArcSP))','UniformOutput',false);
    ArcSP       = [ArcSP(Remain);num2cell(sum(SOIObjX,2))];
end