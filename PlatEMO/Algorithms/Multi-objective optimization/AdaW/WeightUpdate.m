function [Population,W,B] = WeightUpdate(Population,W,Archive,Z,T,Global)
% Weight Update

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Routine to find undeveloped individuals (correspondingly their weights) in the archive set
    % Normalisation
    N_arc         = length(Archive);
    fmin_arc      = min(Archive.objs);
    fmax_arc      = max(Archive.objs);
    Archiveobjs   = (Archive.objs - repmat(fmin_arc,N_arc,1) )./repmat(fmax_arc - fmin_arc,N_arc,1);
    Populaionobjs = (Population.objs - repmat(fmin_arc,Global.N,1) )./repmat(fmax_arc - fmin_arc,Global.N,1);
    % Euclidean distance between individuals in the archive set and individuals in the Population
    dis1 = pdist2(Archiveobjs,Populaionobjs);
    dis1 = sort(dis1,2);
    % Euclidean distance between any two individuals in the archive set
    dis2 = pdist2(Archiveobjs,Archiveobjs);
    dis2 = sort(dis2,2);
    % Calculate the niche size(median of the distances from their closest solution in the archive )
    niche_size = median(dis2(:,2));
    % Find undeveloped 
    Archive_und = Archive(dis1(:,1) >= niche_size);
    N_und = length(Archive_und);
    
    %% If the undeveloped individuals are promising then add them into the evolutionary Population         
    % Obtain their corresponding weights.
	if ~isempty(Archive_und) 
        W1 = (Archive_und.objs - repmat(Z,N_und,1))./repmat( sum(Archive_und.objs,2)-repmat(sum(Z),N_und,1), 1, Global.M );
        for i = 1 : size(W1,1)
            W_all = [W;W1(i,:)];
            B1 = pdist2(W_all,W_all);
            B1(logical(eye(length(B1)))) = inf;
            [~,B1] = sort(B1,2);
            B1 = B1(:,1:T); 

            Population1 = [Population,Archive_und(i)];
            Population2 = Population1(B1(end,:));

            Value_Tche_all = max(abs(Population2.objs-repmat(Z,T,1))./repmat(W1(i,:),T,1),[],2);
            Value_Tche     = max(abs(Archive_und(i).obj -    Z     )./W1(i,:),[],2);
            index = find(Value_Tche_all<Value_Tche, 1);

            if isempty(index)
                % Put the wight into the W, as well as the corresponding solution
                W = [W;W1(i,:)];
                Population = [Population Archive_und(i)];

                % Update neighbour solutions after adding a weight 
                P = B1(end,:);
                g_old = max( abs( Population(P).objs - repmat(Z,T,1) )./W(P,:),[],2 );
                g_new = max( abs( repmat(Archive_und(i).obj,T,1) - repmat(Z,T,1) )./W(P,:),[],2 );
                Population(P(g_old > g_new)) = Archive_und(i);
            end                 
        end
    end
    
    %% Delet the poorly performed weights until the size of W is reduced to N
    % find out the solution that is shared by the most weights in the population
    while length(Population) > Global.N
        [~,ai,bi] = unique(Population.objs,'rows');
        if length(ai) == length(bi)   % If every solution in the population corresponds to only one weight 
            % Normalisation
            fmax  = max(Population.objs,[],1);
            fmin  = min(Population.objs,[],1);
            PCObj = (Population.objs-repmat(fmin,length(Population),1))./repmat(fmax-fmin,length(Population),1);
            % Determine the radius of the niche
            d  = pdist2(PCObj,PCObj);
            d(logical(eye(length(d)))) = inf;
            sd = sort(d,2);
            num_obj = size(Population.objs,2);
            r  = median(sd(:,min(num_obj,size(sd,2))));
            R  = min(d./r,1);
            % Delete solution one by one
            while length(Population) > Global.N
                [~,worst]  = max(1-prod(R,2));
                Population(worst)  = [];
                R(worst,:) = [];
                R(:,worst) = [];
                W(worst,:) = [];
            end
        else
            Index = find(bi==mode(bi));
            Value_Tche2 = max(abs(Population(Index).objs-repmat(Z,size(Index,1),1))./W(Index,:),[],2);
            Index_max= find(Value_Tche2 == max(Value_Tche2));
            Population(Index(Index_max(1)))=[]; 
            W(Index(Index_max(1)),:)=[];             
        end
    end
    % Update the neighbours of each weight
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T); 
end