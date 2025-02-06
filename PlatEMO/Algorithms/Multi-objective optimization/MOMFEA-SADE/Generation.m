function offspring = Generation(Problem,population, Best, DE_Pool, t,RMP,LP,F1,F2,LCR,UCR)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    other_task = 1:length(population);
    other_task(other_task == t) = [];
    for i = 1 : length(population)
        pop_Dec{i} = population{i}.decs;
        pop_Dec{i} = pop_Dec{i}(:,1:end-1);
    end

    x1_task      = other_task(randi(length(other_task)));
    x1_Dec_other = domainAdaption(pop_Dec{t}, pop_Dec{x1_task}, randi(length(population{t}), 1, length(population{t})));
    x2_task      = other_task(randi(length(other_task)));
    x2_Dec_other = domainAdaption(pop_Dec{t}, pop_Dec{x2_task}, randi(length(population{t}), 1, length(population{t})));
    x3_task      = other_task(randi(length(other_task)));
    x3_Dec_other = domainAdaption(pop_Dec{t}, pop_Dec{x3_task}, randi(length(population{t}), 1, length(population{t})));

    x1_Dec_other = [x1_Dec_other,t*ones(size(x1_Dec_other,1),1)];
    x2_Dec_other = [x2_Dec_other,t*ones(size(x2_Dec_other,1),1)];
    x3_Dec_other = [x3_Dec_other,t*ones(size(x3_Dec_other,1),1)];
    for i = 1 : length(population{t})
        offspring(i) = population{t}(i);
        A = randperm(length(population{t}), 4);
        A(A == i) = []; x1 = A(1); x2 = A(2); x3 = A(3);
        CR = LCR + rand() .* (UCR - LCR);

        if rand() < RMP % Random Mating
            x1_Dec = x1_Dec_other(i, :);
            x2_Dec = x2_Dec_other(i, :);
            x3_Dec = x3_Dec_other(i, :);
            switch DE_Pool(i)
                case 1 % DE/best/1/bin
                    offspringDec = population{t}(Best{t}(randi(length(Best{t})))).dec + F1 * (x1_Dec - x2_Dec);
                    offspringDec = DE_Crossover(offspringDec, population{t}(i).dec, CR);
                case 2 % DE/rand/1/bin
                    offspringDec = x1_Dec + F1 * (x2_Dec - x3_Dec);
                    offspringDec = DE_Crossover(offspringDec, population{t}(i).dec, CR);
                case 3 % DE/current-to-rand/1
                    offspringDec = population{t}(i).dec + F2 * (x1_Dec - population{t}(i).dec) + F1 * (x2_Dec - x3_Dec);
            end
        else
            switch DE_Pool(i)
                case 1 % DE/best/1/bin
                    offspringDec = population{t}(Best{t}(randi(length(Best{t})))).dec + F1 * (population{t}(x1).dec - population{t}(x2).dec);
                    offspringDec = DE_Crossover(offspringDec, population{t}(i).dec, CR);
                case 2 % DE/rand/1/bin
                    offspringDec = population{t}(x1).dec + F1 * (population{t}(x2).dec - population{t}(x3).dec);
                    offspringDec = DE_Crossover(offspringDec, population{t}(i).dec, CR);
                case 3 % DE/current-to-rand/1
                    offspringDec = population{t}(i).dec + F2 * (population{t}(x1).dec - population{t}(i).dec) + F1 * (population{t}(x2).dec - population{t}(x3).dec);
            end
        end

        offspringDec(offspringDec > 1) = 1;
        offspringDec(offspringDec < 0) = 0;
        offspringDec(:,end) = t;
        offspring(i) = Problem.Evaluation(offspringDec,[true,DE_Pool(i)]);
    end
end

function OffDec = DE_Crossover(OffDec, ParDec, CR)
    replace = rand(1, size(OffDec, 2)) > CR;
    replace(randi(size(OffDec, 2))) = false;
    OffDec(:, replace) = ParDec(:, replace);
end