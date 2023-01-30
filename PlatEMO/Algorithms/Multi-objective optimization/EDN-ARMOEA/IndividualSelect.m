function xnew = IndividualSelect(PopDec, PopObj, PopMSE, Ke, flag)
% The selection of individuals

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [IDX, ~]=kmeans(PopObj,Ke);
    UIDX=unique(IDX);  
    
    if flag==0
        Uncertainty=mean(PopMSE,2);
        for i=1:length(UIDX)
            Pindex=find(IDX==UIDX(i));
            [~,ind]=max(Uncertainty(Pindex));
            %disp(sprintf('ind is %u', ind)); disp(Pindex');
            xnew(i,:)=PopDec(Pindex(ind),:);          
        end
    else
        Convergence=sqrt(sum(PopObj.^2,2));
        for i=1:length(UIDX)
            Pindex=find(IDX==UIDX(i));
            [~,ind]=min(Convergence(Pindex));
            %disp(sprintf('ind is %u', ind));disp(Pindex');
            xnew(i,:)=PopDec(Pindex(ind),:);      
        end        
    end
end