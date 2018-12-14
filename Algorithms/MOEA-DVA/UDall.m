function Data = UDall(N,M)
% Generate N M-dimensional points by the uniformly design

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is modified from the code in
% http://web.xidian.edu.cn/fliu/lunwen.html

    %% All the possible values of each points
    hm  = find(gcd(1:N,N)==1);
    udt = mod((1:N)'*hm,N); 
    udt(udt==0) = N;
    
    %% Choose M columns among hm as the output, which have the minimum CD2 value
    nCombination = nchoosek(length(hm),M);
    if nCombination < 1e4
        Combination = nchoosek(1:length(hm),M);
        CD2 = zeros(nCombination,1);
        for i = 1 : nCombination
            UT     = udt(:,Combination(i,:));
            CD2(i) = calCD2(UT);
        end
        [~,minIndex] = min(CD2);
        Data = udt(:,Combination(minIndex,:));
    else
        CD2 = zeros(N,1);
        for i = 1 : N
            UT     = mod((1:N)'*i.^(0:M-1),N);
            CD2(i) = calCD2(UT);
        end
        [~,minIndex] = min(CD2);
        Data = mod((1:N)'*minIndex.^(0:M-1),N);
        Data(Data==0) = N;
    end
    Data = (Data-1)/(N-1);
end

function CD2 = calCD2(UT)
% Calculate the CD2 (centered L2-discrepancy) value of the point set, to
% measure the uniformity of the points

    [N,S] = size(UT);
    X     = (2*UT-1)/(2*N);
    
    CS1 = sum(prod(2+abs(X-1/2)-(X-1/2).^2,2));
    CS2 = zeros(N,1);
    for i = 1 : N    
        CS2(i) = sum(prod((1+1/2*abs(repmat(X(i,:),N,1)-1/2)+1/2*abs(X-1/2)-1/2*abs(repmat(X(i,:),N,1)-X)),2));
    end
    CS2 = sum(CS2);
    CD2 = (13/12)^S-2^(1-S)/N*CS1+1/(N^2)*CS2;
end