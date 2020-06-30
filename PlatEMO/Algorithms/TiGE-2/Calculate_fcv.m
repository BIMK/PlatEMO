function fcv = Calculate_fcv(Population)
% calculate normalized  constraints violation(CV) measuring feasibility
    CV_Original=Population.cons;
    CV_Original(CV_Original<=0)=0;
    CV=CV_Original./max(CV_Original);
    CV(:,isnan(CV(1,:))==1)=0;
    fcv=sum(max(0,CV),2)./size(CV_Original,2);
end