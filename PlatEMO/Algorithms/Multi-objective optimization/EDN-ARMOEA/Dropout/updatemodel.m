function net = updatemodel(tr_x, tr_y, Params, net)

    [N,V]=size(tr_x);
    batchsize=V;
    run=8000;
    index=round(1+rand(1,run*batchsize)*(N-1));
    for j=1:run
        x=tr_x(index((j-1)*batchsize+1:j*batchsize),:);
        y=tr_y(index((j-1)*batchsize+1:j*batchsize),:);      
        net=trainNet(x, y, Params, net);
    end
end