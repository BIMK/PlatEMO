function x = NewtonRaphsonMethod(point, precision)
% The Newton-Raphson method

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

    n = size(point, 2);
    x = 1;

    paste_value = x;
    for i=1:100
        % Original function
        f = log(sum(point.^x));

        % Derivative
        numerator = 0;
        denominator = 0;
        for index = 1:n
            if (point(1, index) ~=0)
                numerator = numerator + point(1, index)^x * log(point(1, index));
                denominator = denominator + point(1, index).^x;
            end
        end
        ff = numerator/denominator;

        % zero of function
        x =  x - f /ff;

        if (abs(x-paste_value) <= precision)
            break;
        end
        paste_value = x;
    end
    x = real(x);
end 