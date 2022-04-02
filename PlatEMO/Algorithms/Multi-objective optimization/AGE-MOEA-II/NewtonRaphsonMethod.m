%% Newton-Raphson method
function x = NewtonRaphsonMethod(point, precision)
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