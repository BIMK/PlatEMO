 function Y = sample_dpp(L,k)
% The function of DPPs, the inputs include the kernal matrix L and the size k

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    n     = size(L.D,1);
    [~,i] = sort(L.D);
    v     = i(n-k+1:n);
    k     = length(v);
    V     = L.V(:,v);

    % iterate
    Y = zeros(k,1);
    for i = k:-1:1

        if size(V,2)==1 && i~=1
            V=abs(V);
            [~,b1]=sort(V);
            bb=b1(n-i+1:n);
            for ii=i:-1:1
                Y(ii)=bb(ii);
            end
            break;

        end
        % compute probabilities for each item
        P = sum(V.^2,2);
        P = P / sum(P);

        [~,Y(i)] = max(P);

        % choose a vector to eliminate
        j  = find(V(Y(i),:),1);
        Vj = V(:,j);
        V  = V(:,[1:j-1 j+1:end]);

        % update V

        V = V - bsxfun(@times,Vj,V(Y(i),:)/Vj(Y(i)));

        % orthogonalize
        V = orth(V);
    end
    Y = sort(Y);
end