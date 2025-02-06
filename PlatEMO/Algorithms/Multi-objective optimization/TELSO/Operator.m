function [Offspring,Mask] = Operator(Population,Mask,Problem)
% The learning swarm optimizer of TELSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Sheng Qi (email: 2745679162@qq.com)

    %% Parameter setting
    PopulationDec = Population.decs;
    [N,D]         = size(PopulationDec);
	PopulationVel = Population.adds(zeros(N,D));
    OffDec        = [];
    OffVel        = []; 
    
    %% Refactoring the real variable    
    for i = 1 : (N-2)
        rand_indices = randi([i+1,N],1,2);
        selected_particles = Population(rand_indices); 
        r1       = repmat(rand(1,1),1,D);
        r2       = repmat(rand(1,1),1,D);
        r3       = repmat(rand(1,1),1,D);
        p_OffVel = r1.*PopulationVel(i)+r2.*(selected_particles(1).dec-Population(i).dec)+r3.*(selected_particles(2).dec-Population(i).dec);
        p_OffDec = Population(i).dec+p_OffVel;
        OffDec   = [OffDec;p_OffDec];
        OffVel   = [OffVel;p_OffVel];

        p1_mask = Mask(rand_indices(1), :);
        p2_mask = Mask(rand_indices(2), :);
        if rand < 0.5
            same_elems = (p1_mask == p2_mask & p2_mask ==1);
        else
             same_elems = (p1_mask == p2_mask & p2_mask ==0);
        end
        same_cols = find(same_elems);
        D_number  = numel(same_cols);  
        keep      = floor(D_number * Problem.FE/Problem.maxFE);  
        same_cols = same_cols(1:keep);  
        if same_cols ~= 0
            for j = same_cols
                Mask(i,j) = Mask(rand_indices(1),j); 
            end
        end
    end

	%% Add the two winners
    OffDec    = [OffDec;Population(N-1).dec;Population(N).dec];
    OffVel    = [OffVel;Population(N-1).add;Population(N).add];
    Offspring = Problem.Evaluation(OffDec.*Mask,OffVel);
end