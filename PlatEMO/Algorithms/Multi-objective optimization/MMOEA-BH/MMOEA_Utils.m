classdef MMOEA_Utils
% MMOEA_Utils - Static class for all utility and auxiliary functions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    methods(Static)
        %% 1. APC Clustering
        function varargout = APC(data)
            N    = size(data,1);
            DM   = pdist2(data,data);
            S    = -DM;
            meds = median(S,'all');
            n    = 5;
            
            S(logical(eye(N))) = n*meds;
            A     = zeros(N,N);
            R     = zeros(N,N); 
            lam   = 0.6; 
            count = 0;
            
            for iter = 1 : 500
                Eold  = A + R;
                Rold  = R;
                AS    = A+S;
                [Y,I] = max(AS,[],2);
                for i = 1 : N
                    AS(i,I(i)) = -inf;
                end
                Y2 = max(AS,[],2);
                R  = S-Y;
                for i = 1 : N
                    R(i,I(i)) = S(i,I(i))-Y2(i);
                end
                R = (1-lam)*R+lam*Rold; 
                
                Aold = A;
                Rp   = max(R,0);    
                Rp(logical(eye(N))) = R(logical(eye(N)));
                A    = sum(Rp,1)-Rp;
                dA   = diag(A);
                A    = min(A,0);    
                A(logical(eye(N))) = dA;
                A    = (1-lam)*A+lam*Aold; 
                E    = A + R;
                
                if diag(Eold) == diag(E)
                    count = count + 1;
                    if count == 10
                        break;
                    end
                else
                    count=0;
                end
            end
            E      = R + A; 
            I      = find(diag(E)>0);
            K      = length(I);
            [~, c] = max(E(:,I),[],2);
            c(I)   = 1 : K;                 
            idx    = I(c); 
            
            maxCluster = 0;
            for i = unique(idx)'
                maxCluster  = maxCluster + 1;
                idx(idx==i) = maxCluster;
            end
            varargout = {idx, maxCluster};
        end
        
        %% 2. Calculate PCCS
        function varargout = calc_pccs(obj)
            K = size(obj,1);
            if K >= 2
                fmax    = max(obj);
                fmin    = min(obj);
                L       = ceil(K * (obj-fmin)./(fmax-fmin));
                L(L==0) = L(L==0) + 1;
                
                PCD = pdist2(L,L,'cityblock');
                PCD(logical(eye(K))) = [];
                PCD         = reshape(PCD,K-1,K)';
                PCD(PCD==0) = PCD(PCD==0) + 0.5;
                density     = sum(1./PCD.^2,2);
                varargout   = {density};
            else
                varargout   = {1};
            end
        end
        
        %% 3. Non-dominated PCCS Sort
        function [selected_population, sorted_population, sorted_SCD]= nd_pccs_sort(population, aSize)
            if nargin < 2
                aSize = 1000;
            end
            PopDec = population.decs;
            PopObj = population.objs;
            [N,~]  = size(PopDec);
            if N == 1
                selected_population = population;
                sorted_population   = population;
                sorted_SCD          = 1;
            else        
                CD_x    = normalize( -MMOEA_Utils.calc_pccs(PopDec), 'range' );
                CD_f    = normalize( -MMOEA_Utils.calc_pccs(PopObj), 'range' );
                idx_max = bitor(CD_x>mean(CD_x), CD_f>mean(CD_f));
                idx_min = ~idx_max;
                
                CD(idx_max) = CD_x(idx_max) + CD_f(idx_max);
                CD(idx_min) = min(CD_x(idx_min), CD_f(idx_min));
                [sorted_SCD, idx]   = sort(CD,'descend');
                selected_population = population(idx(1));
                
                if numel(idx) > aSize
                    sorted_population = population(idx(1:aSize));
                    sorted_SCD        = sorted_SCD(1:aSize);
                else
                    sorted_population = population(idx);
                end
            end
        end
        
        %% 4. REP Selection for PSO
        function [pbest, lbest] = rep_selection_pso(Problem, PBA, PBA_SCD, LBA, LBA_SCD, idx_APC, maxCluster_APC)
            for i = 1 : numel(PBA)
                idx1     = 1;
                pbest(i) = PBA{i}(idx1);
            end
            if numel(idx_APC)<Problem.N
                k       = ceil(Problem.N/numel(idx_APC));
                idx_APC = repmat(idx_APC,k,1);    
            end
            idx_APC = idx_APC(1:Problem.N);
            for i = 1 : maxCluster_APC
                lbest_temp = LBA{i};
                pos        = i==idx_APC;
                idx        = TournamentSelection(2,sum(pos),-LBA_SCD{i});
                lbest(pos) = lbest_temp(idx);
            end
        end
        
        %% 5. Get Original Volume
        function [orig_vol, center] = get_orig_vol(Problem, maxCluster, LBA)
            for i = 1 : maxCluster
                pop         = LBA{i};
                pop_dec     = pop.decs;
                center(i,:) = mean(pop_dec,1);
            end
            grid_num = ceil(power(maxCluster,1/Problem.D));
            len1     = (Problem.upper-Problem.lower)/grid_num;
            orig_vol = prod(len1);
        end
        
        %% 6. Get Upper and Lower bounds
        function [new_up_s, new_low_s] = get_upper_lower(Problem, pop, orig_vol)
            pop_dec          = pop.decs;
            upper_small      = max(pop_dec,[],1);
            lower_small      = min(pop_dec,[],1);
            dimensions_small = upper_small - lower_small;
            
            ratio          = dimensions_small ./ max(dimensions_small);
            adjusted_ratio = max(ratio, 0.2);
            ratio          = adjusted_ratio ./ sum(adjusted_ratio);
            
            new_volume     = orig_vol;
            new_dimensions = (new_volume / prod(ratio))^(1/length(ratio)) * ratio;
            
            center    = (upper_small + lower_small) / 2;
            new_low_s = center - new_dimensions / 2;
            new_up_s  = new_low_s + new_dimensions;
            new_up_s  = max(min(new_up_s,Problem.upper),Problem.lower);
            new_low_s = max(min(new_low_s,Problem.upper),Problem.lower);
        end
        
        %% 7. Clip Upper and Lower bounds
        function [new_new_up_s, new_new_low_s] = clip_upper_lower(new_up_s, new_low_s, center)
            inside_each_dim  = (center >= new_low_s) & (center <= new_up_s);
            inside_rectangle = all(inside_each_dim, 2);
            points_inside    = center(inside_rectangle, :);
            num_inside       = sum(inside_rectangle); 
            
            if num_inside < 2
                new_new_up_s  = new_up_s;
                new_new_low_s = new_low_s;
            else
                original_volume = prod(new_up_s - new_low_s);
                target_volume   = original_volume / num_inside;
                X               = points_inside(:, 1:end-1);
                y               = points_inside(:, end);
                mdl             = fitlm(X, y);
                
                coefficients        = round(1./abs([mdl.Coefficients.Estimate(2:end)', 1]), 2);
                coefficients(isinf(coefficients)) = 10^6;
                original_dimensions = new_up_s - new_low_s;
                syms a
                equation_terms = 1;
                for i = 1 : length(original_dimensions)
                    equation_terms = equation_terms * (original_dimensions(i) - coefficients(i) * a);
                end
                
                equation  = equation_terms == target_volume;
                solutions = double(solve(equation, a));
                valid_solution = [];
                for i = 1 : length(solutions)
                    if all(solutions(i) > 0) && all(original_dimensions ./ coefficients > solutions(i))
                        valid_solution = solutions(i);
                        break;
                    end
                end
                
                if isempty(valid_solution)
                    error('No valid solution found');
                end
                adjusted_target_side_length = original_dimensions - coefficients * valid_solution;
                centroid      = mean((new_up_s + new_low_s) / 2, 1);
                new_new_low_s = centroid - adjusted_target_side_length / 2;
                new_new_up_s  = centroid + adjusted_target_side_length / 2;
                new_new_up_s  = max(new_new_up_s, new_low_s);
                new_new_low_s = min(new_new_low_s, new_up_s);
            end
        end
    end
end