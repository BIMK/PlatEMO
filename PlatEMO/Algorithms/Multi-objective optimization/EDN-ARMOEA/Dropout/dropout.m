function [x, index1, index2] = dropout(x,dropP)

    [n1, n2]=size(x);
    array=rand(n1, n2)<repmat(dropP, n1, n2);
    x(find(array==1))=0;
    [index1, index2]=find(array==0);
    %[index3, index4]=find(array==1);
    x = x/(1-dropP);
end