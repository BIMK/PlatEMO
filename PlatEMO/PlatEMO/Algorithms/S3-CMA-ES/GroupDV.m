function group = GroupDV(Global,DV,PV,nPerGroup)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen

   ub  = Global.upper(DV);
   lb  = Global.lower(DV); 
   dim = length(DV);	% the dim of the distance variables
   
   fixPV  = (Global.lower(PV)+Global.upper(PV))/2;	% the value vector for position variables
   meanDV = (ub+lb)/2;
   
   f_archiveF = repmat(zeros(dim,dim),1,Global.M);	% cell(dim, dim);
   
   fhatDVDecs     = repmat(lb,dim,1) + eye(dim).*repmat(meanDV,dim,1);
   fhatDecs       = [repmat(fixPV,dim,1) fhatDVDecs];
   fhatPopulation = INDIVIDUAL(fhatDecs);
   fhat_archiveF  = fhatPopulation.objs;
   
   lambdaF = repmat(zeros(dim,dim),1,Global.M);	% cell(dim, dim);
   
   p1Dec   = [fixPV lb]; 
   tempFp1 = INDIVIDUAL(p1Dec);
   fp1F    = tempFp1.objs;
   
   for i = 1 : dim-1      
       fp2 = fhat_archiveF(i,:);
       for j = i+1 : dim
           
           fp3 = fhat_archiveF(j,:);

           p4      = lb;
           p4(i)   = meanDV(i);	% temp;
           p4(j)   = meanDV(j);	% temp;
           p4Dec   = [fixPV p4];
           tempFp4 = INDIVIDUAL(p4Dec);
           fp4     = tempFp4.objs;
           
           f_archiveF(i,j:dim:size(f_archiveF,2)) = fp4;

           d1 = fp2 - fp1F;
           d2 = fp4 - fp3;
           
           lambdaF(i,j:dim:size(lambdaF,2)) = abs(d1-d2);
       end
   end
   
   % Check for each objective
   bigTheta = false(length(DV));
   
   for i = 1 : Global.M
       
       lambda = lambdaF(:,(i-1)*dim+1:i*dim);
       fp1    = fp1F(i);
       
       fhat_archive = fhat_archiveF(:,i);
       
       tempF_archive = f_archiveF(:,(i-1)*dim+1:i*dim);
       f_archive     = tempF_archive + tempF_archive';
       
       F1 = ones(dim,dim)*fp1;
       F2 = repmat(fhat_archive',dim,1);
       F3 = repmat(fhat_archive,1,dim);
       F4 = f_archive;
       
       FS   = cat(3,F1,F2,F3,F4);
       Fmax = max(FS,[],3);
       
       FS       = cat(3,F1+F4,F2+F3);
       Fmax_inf = max(FS,[],3);

       theta = false(dim);

       muM       = eps/2;
       gamma     = @(n)((n.*muM)./(1-n.*muM));
       errlb     = gamma(2)*Fmax_inf;
       errub     = gamma(dim^0.5)*Fmax;
       I2        = lambda >= errub;
       theta(I2) = 1;
       
       % add, then find not less than 1
       bigTheta = bigTheta + theta;
   end
   
   bigTheta(bigTheta>0) = true;

   L = size(bigTheta,1);	% number of vertex
   
   % Breadth-first search
   labels = zeros(1,L);     % all vertex unexplored at the begining
   rts    = [];
   ccc    = 0;              % connected components counter
   while true
       ind = find(labels==0);
       if ~isempty(ind)
           fue  = ind(1);   % first unexplored vertex
           rts  = [rts fue];
           list = [fue];
           ccc  = ccc + 1;
           labels(fue) = ccc;
           while true
               list_new = [];
               for lc = 1 : length(list)
                   p   = list(lc);              % point
                   cp  = find(bigTheta(p,:));	% points connected to p
                   cp1 = cp(labels(cp)==0);     % get only unexplored vertecies
                   labels(cp1) = ccc;
                   list_new    = [list_new cp1];
               end
               list = list_new;
               if isempty(list)
                   break;
               end
           end
       else
           break;
       end
   end
   
   group_num = max(labels);
   allgroups = cell(1,group_num);
   for i = 1 : group_num
       allgroups{i} = find(labels==i);
   end
   
   h = @(x)(length(x)==1);
   sizeone = cellfun(h,allgroups);
   
   seps = allgroups(sizeone);
   seps = cell2mat(seps);
   
   allgroups(sizeone) = [];
   nonseps            = allgroups;

   group = {};
   for ns = 1 : length(nonseps)     % the non-seperate variables
       group = {group{1:end} DV(nonseps{ns})};
   end

   for g = 1 : nPerGroup : length(seps)
       index = seps(g:min(g+nPerGroup-1,length(seps)));
       group = {group{1: end} DV(index)};
   end
end