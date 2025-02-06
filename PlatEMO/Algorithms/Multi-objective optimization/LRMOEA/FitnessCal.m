function[score] = FitnessCal(Mask,FrontNo,score)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,d] = size(Mask);
    Fit   = zeros(1,d);
    b     = 0;
    h     = 0;
    for i = 1 : N
        h = h + 1;
        if FrontNo(h) <= 1
            b = b + 1;
            a = 1;
            for j = 1 : d
                fit = 1/(2+sqrt(N)*Mask(i,j));
                Fit(1,a)  = fit ;
                a = a + 1;
            end
            score = score + Fit;
        end
    end
    score = score./b;
end