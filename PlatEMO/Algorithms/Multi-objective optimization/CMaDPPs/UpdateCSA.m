function [CSA]= UpdateCSA(CSA,Offspring,MaxSize,epsilon)
% Update CSA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    % Select epsilon feasible solutions
    CSA = [Offspring,CSA];
    CV = sum(max(0,CSA.cons),2);
    fIndex = all(CV <= epsilon,2);
    CSA = CSA(fIndex);
	CSAObj = CSA.objs;
    
    [N,M] = size(CSAObj);

    if N <= MaxSize
        return;
    end
    
    cunum=ceil((MaxSize)/(3*M));
    DAObj2 = CSAObj;
    DAObj3=DAObj2.^2;
    DAO=sum(DAObj3,2);
    CHIndex = [];
    minIndex2=[];
    for i=1:M
        minfx=CSAObj(:,i);
        tempminfx=min(minfx)+1E-6;

        if sum(minfx<tempminfx)>cunum
            mindex=find((minfx<tempminfx));
            [~,minIndex]=sort(DAO(mindex));
            minIndex=mindex(minIndex(1:cunum));
        else
            [~,minIndex]=sort(minfx);
            minIndex=minIndex(1:cunum);
        end
      
        minfx2=DAO-DAObj3(:,i);
        tempminfx=min(minfx2)+1E-6;

        if sum(minfx2<tempminfx)>2*cunum
            mindex=find((minfx2<tempminfx));
            [~,minndex]=sort(DAO(mindex));
            index=mindex(minndex(1:2*cunum));
        else
            [~,index]=sort(minfx2);
            index=index(1:2*cunum);
        end

       CHIndex = [CHIndex,index'];

       minIndex2=[minIndex2,minIndex'];
    end

    Choose = [CHIndex,minIndex2];
    CSA    = CSA(Choose);
end