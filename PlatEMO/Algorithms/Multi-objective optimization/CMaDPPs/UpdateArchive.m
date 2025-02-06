function Population = UpdateArchive(PopOff,MaxSize,z,znad,DAobj,theta,epsilon)
% Update the archive for feasible solutions by MaOEADPPs approach

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    PopObjs = PopOff.objs;
    PopCons = PopOff.cons;
    CV = sum(max(0,PopCons),2);
    fIndex = all(CV <= epsilon,2);
    PopCons(fIndex,:) = 0;
    
    ND = NDSort(PopObjs,PopCons,1);
    Population = PopOff(ND==1);
    
    CAObj=Population.objs;
    [~,ia,~] = unique(CAObj,'rows');
    Population = Population(ia);

    N  = length(Population);
    
    if N <= MaxSize
        return;
    end
    
    CAObj=Population.objs;
    CAObj2 = (CAObj-repmat(z,N,1))./(repmat(znad,N,1)-repmat(z,N,1));
    
    nad=max(DAobj)+1E-6;
    nn=sum((nad-CAObj)<0,2);
    [N1,~]=size(DAobj);
    DAobj = (DAobj-repmat(z,N1,1))./(repmat(znad,N1,1)-repmat(z,N1,1));

    big=max(sqrt(sum(DAobj.^2,2)))+1E-6;
    nbig=sqrt(sum(CAObj2.^2,2))>big;
    nn=nn+nbig;
   
    D = pdist2(CAObj2,CAObj2,'cosine');

    D(1:size(D,1)+1:end) = 0;

    H=(1-theta)*exp(-D);

    value=(sum((CAObj2).^2,2));

    value(nn==0)=min(value)/2;
    value=value/max(value);
    L=value*value';
    H=H./L;
    
    L = decompose_kernel(H);
    Choose = sample_dpp(L,MaxSize);
    Population = Population(Choose);
end