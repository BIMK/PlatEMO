function f = Distance_function(x, problem)
warning off
upper = repmat(10,size(x,1),size(x,2));
lower = repmat(0,size(x,1),size(x,2));
x = (upper - lower).* x + lower;

if problem  ==1
    % Sphere function
    f = sum(x.^2,2);
elseif problem == 2
    % Schwefel function
     f = max(abs(x),[],2);
elseif problem == 3
    %  Rastrigin function
    f = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
elseif problem == 4
    % Griewank function
    f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
elseif problem == 5
    % Ackley function
     f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
end
f(f<1e-8) = 0;
end

