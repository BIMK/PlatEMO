function [PopObj, PopStd] = Estimate(PopDec, net, Params, M)

    %tau=Params.tau;
    %interval=Params.interval;
    % Calculate the objective values according to the decision
    % variables, note that here the objective values of multiple
    % solutions are calculated at the same time
    ps=Params.ps;qs=Params.qs;
    x=mapminmax('apply',PopDec',ps);x=x';


        for i=1:100
            sum_y=testNet(x, net, Params);
            sum_y=mapminmax('reverse',sum_y',qs);sum_y=sum_y';
            sum_ysq=sum_y.^2;
            result(i).y=sum_y;
            result(i).ysq=sum_ysq;
        end

    %Take the result of the most recent interval
    array1=[result.y];array2=[result.ysq];
    for i=1:M
        mu(:,i)=mean(array1(:,i:M:end),2);
        s2(:,i)=mean(array2(:,i:M:end),2);
    end
    std=sqrt(s2-mu.^2);

    %alpha=2;
    PopObj=mu;%-alpha*std;
    PopStd=std;
end