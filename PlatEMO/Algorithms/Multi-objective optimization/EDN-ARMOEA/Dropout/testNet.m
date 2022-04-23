function x7 = testNet(x, net, Params)

    N=size(x,1);
    W=net.W;B=net.B;
    dropP=Params.dropP;
    %% forward
    [x1,~]=dropout(x,dropP(1));%N*V
    x2=x1*W{1}+repmat(B{1},N,1);%N*neuronN
    x3=max(0,x2);%ReLU
    [x4, ~]=dropout(x3,dropP(2));
    x5=x4*W{2}+repmat(B{2},N,1);
    x6=(exp(x5)-exp(-x5))./(exp(x5)+exp(-x5));%Sigmoid%N*neuronN
    x7=x6*W{3}+repmat(B{3},N,1);%N*M
end