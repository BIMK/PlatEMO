function net = trainNet(x, y, Params, net)

    N=size(x,1);
    dropP=Params.dropP;learnR=Params.learnR;
    W=net.W;B=net.B;

    %% forward
    [x1,~]=dropout(x,dropP(1));%N*V
    x2=x1*W{1}+repmat(B{1},N,1);%N*neuronN
    x3=max(0,x2);%ReLU
    [x4, index1, index2]=dropout(x3,dropP(2));
    x5=x4*W{2}+repmat(B{2},N,1);
    x6=(exp(x5)-exp(-x5))./(exp(x5)+exp(-x5));%Sigmoid%N*neuronN
    x7=x6*W{3}+repmat(B{3},N,1);%N*M

    %% backward
    %cost_loss=sum(sum(0.5*(x7-y).^2));
    %fc x6-x7
    e=x7-y;%N*M
    dW{3}=x6'*e;%neuronN*M
    dB{3}=sum(e);%1*M
    dx6=e*W{3}';%N*neuronN
    %sigmoid x5-x6
    dx5=dx6.*(1-x6.^2);%N*neuronN
    %fc x4-x5
    dW{2}=x4'*dx5;%neuronN*neuronN
    dB{2}=sum(dx5);%1*neuronN
    dx4=dx5*W{2}';%N*neuronN
    %dp x3-x4
    %tic;
    dx3=zeros(size(x3));%N*neuronN
    for i=1:length(index1)
        dx3(index1(i),index2(i))=dx4(index1(i), index2(i))/(1-dropP(2));
    end
    %ReLU x2-x3
    dx2=dx3;%N*neuronN
    dx2(find(x3<=0))=0;%N*neuronN
    %fc x1-x2
    dW{1}=x1'*dx2;%V*neuronN
    dB{1}=sum(dx2);%1*neuronN

    decay=1e-05;
    for k=1:3
        W{k}=W{k}-(decay*W{k}+dW{k})/N*learnR;
        B{k}=B{k}-dB{k}/N*learnR;
    end
    net.W=W;net.B=B;
end