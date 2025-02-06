function offspring = Generation(Problem,population, rank, model, t,Tau1,Tau2)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    for i = 1 : length(population{t})
        offspring(i) = population{t}(i);
        
        % Parameter disturbance
        F  = normrnd(population{t}(i).add(1), 0.1);
        F  = min(max(F, 0.2), 1.2);
        CR = normrnd(population{t}(i).add(2), 0.1);
        CR = min(max(CR, 0), 1);
        TR = normrnd(population{t}(i).add(3), 0.1);
        TR = min(max(TR, 0), 1);
        KP = normrnd(population{t}(i).add(4), 0.1);
        if KP < 0
            KP = 1 + KP;
        elseif KP > 1
            KP = KP - 1;
        end
    
        % Parameter mutation
        if rand() < Tau1
            F = 0.2 + rand();
        end
        if rand() < Tau1
            CR = rand();
        end
        if rand() < Tau2
            TR = rand();
        end
        if rand() < Tau2
            KP = rand();
        end
    
        % Select individuals (rank-DE)
        Np = length(population{t});
        x1 = randi(Np);
        while rand() > (Np - rank{t}(x1)) / Np || x1 == i
            x1 = randi(Np);
        end
        x2 = randi(Np);
        while rand() > (Np - rank{t}(x2)) / Np || x2 == i || x2 == x1
            x2 = randi(Np);
        end
        x3 = randi(Np);
        while x3 == i || x3 == x1 || x3 == x2
            x3 = randi(Np);
        end
        xDeci = population{t}(i).dec;
        xDec1 = population{t}(x1).dec;
        xDec2 = population{t}(x2).dec;
        xDec3 = population{t}(x3).dec;
    
        % Knowledge transfer
        if rand() < TR
            k = randi(length(population)); % help task
            while (k == t), k = randi(length(population)); end
            Np = length(population{k});
    
            rnd = KP;
            if rnd > 1/2
                xDeck = population{k}(randi(Np)).dec;
            else
                xDeck = population{k}(randi(Np)).dec;
                xDeck = (xDeck -model{k}.mean) ./ model{k}.std;
                xDeck = model{t}.mean + model{t}.std .* xDeck;
            end
            xDec2 = xDeck;
        end
    
        offspringdec = xDec1 + F * (xDec2 - xDec3);
        offspringdec = DE_Crossover(offspringdec, xDeci, CR);
        offspringdec = min(max(offspringdec, 0), 1);
        offspringdec(end) = t;
        offspring(i) = Problem.Evaluation(offspringdec,[F,CR,TR,KP]);
    end
end

function OffDec = DE_Crossover(OffDec, ParDec, CR)
    replace = rand(1, size(OffDec, 2)) > CR;
    replace(randi(size(OffDec, 2))) = false;
    OffDec(:, replace) = ParDec(:, replace);
end