function Offspring1=DEgenerator(Population1,Population2,Problem,Fitness1)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    p1 = Population1.decs;
    lu = [Problem.lower;Problem.upper];
    [popsize1,n] = size(p1);
    trial = zeros(popsize1,n);
    p2    = Population2.decs;
    
    for i = 1 : popsize1
        l = rand;
        if l <= 1/3
            F = 0.6;
        elseif l <= 2/3
            F = 0.8;
        else
            F = 1.0;
        end
    
        l = rand;
        if l <= 1/3
            CR = 0.1;
        elseif l <= 2/3
            CR = 0.2;
        else
            CR = 1.0;
        end
        indexset    = 1 : popsize1;
        indexset(i) = [];
        r1  = floor(rand*(popsize1-1))+1;
        xr1 = indexset(r1);
        indexset(r1) = [];
        r2  = floor(rand*(popsize1-2))+1;
        xr2 = indexset(r2);
        indexset(r2) = [];
        r3  = floor(rand*(popsize1-3))+1;
        xr3 = indexset(r3);
    
        Fr    = find(Fitness1 < 1);
        rr    = floor(rand*length(Fr))+1;
        best1 = Fr(rr);
    
        o = rand;
    
        if o < 1/2
            v = p1(i,:) + F*(p2(xr1,:)-p1(i,:)) + F*(p1(xr2,:)-p1(xr3,:));
        else
            v = p1(i,:) + F*(p1(best1,:)-p1(i,:)) + F*(p1(xr2,:)-p1(xr3,:));
        end
        % Handle the elements of the mutant vector which violate the boundary
        w = find(v < lu(1,:));
        if ~isempty(w)
            v(1, w) = 2 * lu(1, w) -  v(1, w);
            w1 = find( v(1, w) > lu(2, w));
            if ~isempty(w1)
                v(1, w(w1)) = lu(2, w(w1));
            end
        end
        y = find(v > lu(2, :));
        if ~isempty(y)
            v(1, y) =  2 * lu(2, y) - v(1, y);
            y1 = find(v(1, y) < lu(1, y));
            if ~isempty(y1)
                v(1, y(y1)) = lu(1, y(y1));
            end
        end
        % Binomial crossover
        t = rand(1, n) < CR;
        j_rand = floor(rand * n) + 1;
        t(1, j_rand) = 1;
        t_ = 1 - t;
        trial(i, :) = t .* v + t_ .* p1(i,:);
    end
    Offspring1 = Problem.Evaluation(trial);
end