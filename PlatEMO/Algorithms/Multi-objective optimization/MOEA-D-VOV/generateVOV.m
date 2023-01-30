function VOV = generateVOV(obj,theta)
% Generate the virtual objective vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    obj    = normalize(obj,'range');
    [N,M]  = size(obj);
    [W,N2] = UniformPoint(2e4,M,'ILD');
    VOV    = zeros(N2,M);
    flag   = false(N2,1);

    normW   = sqrt(sum(W.^2,2));
    normObj = sqrt(sum(obj.^2,2));
    for i = 1 : N2
        CosineVOV = sum(obj.*repmat(W(i,:),N,1),2)./normW(i,:)./normObj;
        d2 = normObj.*sqrt(1-CosineVOV.^2);
        [mind2,I] = min(d2);
        d1 = normObj(I)*CosineVOV(I);
        r = d1/norm(W(i,:));

        VOV(i,:) = W(i,:).*r;
        flag(i)  = mind2 < theta;
    end
    VOV = VOV(flag,:);
end 