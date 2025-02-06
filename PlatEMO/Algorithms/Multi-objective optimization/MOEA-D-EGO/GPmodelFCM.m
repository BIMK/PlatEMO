function [model,centers] = GPmodelFCM(train_x,train_y,L1,L2)
% Fuzzy clustering-based method for modeling c_size* M models, where c_size
% is the number of clusters and M the number of objectives.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao (liazhao5-c@my.cityu.edu.hk).

    [K,M] = size(train_y);
    D     = size(train_x,2);
    if K <= L1
        % all the K evaluated points are directly used for building GP model
        csize = 1;
        [centers,~] = FCM(train_x,csize,[2 NaN 0.05 false]);
        model = cell(1,M);
        theta = (K ^ (-1 ./ K)) .* ones(1, D);
        for j = 1 : M
            model{1,j}= Dacefit(train_x,train_y(:,j),'regpoly0','corrgauss',theta,1e-6*ones(1,D),20*ones(1,D));
        end 
    else
        % FuzzyCM
        csize       = 1 + ceil((K-L1)/L2);
        [centers,~] = FCM(train_x,csize,[2 NaN 0.05 false]);
        dis         = pdist2(train_x,centers);
        [~,index]   = sort(dis);
        theta       = (L1 ^ (-1 ./ L1)) .* ones(1, D);

        %% Build GP model of each objective for each cluster 
        model   = cell(csize, M);
        for i   = 1 : csize
            for j = 1 :  M
                temp_index = index(1:L1,i);
                model{i,j} = Dacefit(train_x(temp_index,:),train_y(temp_index,j), 'regpoly0','corrgauss', theta, 1e-6.*ones(1,D), 20.*ones(1,D));
            end
        end
    end

end
% >>>>>>>>>>>>>>>>   Auxiliary functions  ==================== 
function [center, U, obj_fcn] = FCM(data, cluster_n, options)
    %FCM Data set clustering using fuzzy c-means clustering.
    %
    %   [CENTER, U, OBJ_FCN] = FCM(DATA, N_CLUSTER) finds N_CLUSTER number of
    %   clusters in the data set DATA. DATA is size M-by-N, where M is the number of
    %   data points and N is the number of coordinates for each data point. The
    %   coordinates for each cluster center are returned in the rows of the matrix
    %   CENTER. The membership function matrix U contains the grade of membership of
    %   each DATA point in each cluster. The values 0 and 1 indicate no membership
    %   and full membership respectively. Grades between 0 and 1 indicate that the
    %   data point has partial membership in a cluster. At each iteration, an
    %   objective function is minimized to find the best location for the clusters
    %   and its values are returned in OBJ_FCN.
    %
    %   [CENTER, ...] = FCM(DATA,N_CLUSTER,OPTIONS) specifies a vector of options
    %   for the clustering process:
    %       OPTIONS(1): exponent for the matrix U             (default: 2.0)
    %       OPTIONS(2): maximum number of iterations          (default: 100)
    %       OPTIONS(3): minimum amount of improvement         (default: 1e-5)
    %       OPTIONS(4): info display during iteration         (default: 1)
    %   The clustering process stops when the maximum number of iterations
    %   is reached, or when the objective function improvement between two
    %   consecutive iterations is less than the minimum amount of improvement
    %   specified. Use NaN to select the default value.
    %
    %   Example
    %       data = rand(100,2);
    %       [center,U,obj_fcn] = fcm(data,2);
    %       plot(data(:,1), data(:,2),'o');
    %       hold on;
    %       maxU = max(U);
    %       % Find the data points with highest grade of membership in cluster 1
    %       index1 = find(U(1,:) == maxU);
    %       % Find the data points with highest grade of membership in cluster 2
    %       index2 = find(U(2,:) == maxU);
    %       line(data(index1,1),data(index1,2),'marker','*','color','g');
    %       line(data(index2,1),data(index2,2),'marker','*','color','r');
    %       % Plot the cluster centers
    %       plot([center([1 2],1)],[center([1 2],2)],'*','color','k')
    %       hold off;
    %
    %   See also FCMDEMO, INITFCM, IRISFCM, DISTFCM, STEPFCM.
    
    %   Roger Jang, 12-13-94, N. Hickey 04-16-01
    %   Copyright 1994-2018 The MathWorks, Inc. 
    
    if nargin ~= 2 && nargin ~= 3
	    error(message("fuzzy:general:errFLT_incorrectNumInputArguments"))
    end
    
    data_n = size(data, 1);
    
    % Change the following to set default options
    default_options = [2;	% exponent for the partition matrix U
		    100;	% max. number of iteration
		    1e-5;	% min. amount of improvement
		    1];	% info display during iteration 
    
    if nargin == 2
	    options = default_options;
    else
	    % If "options" is not fully specified, pad it with default values.
	    if length(options) < 4
		    tmp = default_options;
		    tmp(1:length(options)) = options;
		    options = tmp;
	    end
	    % If some entries of "options" are nan's, replace them with defaults.
	    nan_index = find(isnan(options)==1);
	    options(nan_index) = default_options(nan_index);
	    if options(1) <= 1
		    error(message("fuzzy:general:errFcm_expMustBeGtOne"))
	    end
    end
    
    expo = options(1);		% Exponent for U
    max_iter = options(2);		% Max. iteration
    min_impro = options(3);		% Min. improvement
    display = options(4);		% Display info or not
    
    obj_fcn = zeros(max_iter, 1);	% Array for objective function
    
    U = initfcm(cluster_n, data_n);			% Initial fuzzy partition
    % Main loop
    for i = 1:max_iter
	    [U, center, obj_fcn(i)] = stepfcm(data, U, cluster_n, expo);
	    if display
		    fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	    end
	    % check termination condition
	    if i > 1
		    if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end
	    end
    end
    
    iter_n = i;	% Actual number of iterations 
    obj_fcn(iter_n+1:max_iter) = [];
end

function out = distfcm(center, data)
    %DISTFCM Distance measure in fuzzy c-mean clustering.
    %	OUT = DISTFCM(CENTER, DATA) calculates the Euclidean distance
    %	between each row in CENTER and each row in DATA, and returns a
    %	distance matrix OUT of size M by N, where M and N are row
    %	dimensions of CENTER and DATA, respectively, and OUT(I, J) is
    %	the distance between CENTER(I,:) and DATA(J,:).
    %
    %       See also FCMDEMO, INITFCM, IRISFCM, STEPFCM, and FCM.
    
    %	Roger Jang, 11-22-94, 6-27-95.
    %       Copyright 1994-2016 The MathWorks, Inc. 
    
    out = zeros(size(center, 1), size(data, 1));
    
    % fill the output matrix
    
    if size(center, 2) > 1
        for k = 1:size(center, 1)
	    out(k, :) = sqrt(sum(((data-ones(size(data, 1), 1)*center(k, :)).^2), 2));
        end
    else	% 1-D data
        for k = 1:size(center, 1)
	    out(k, :) = abs(center(k)-data)';
        end
    end
end
function U = initfcm(cluster_n, data_n)
    %INITFCM Generate initial fuzzy partition matrix for fuzzy c-means clustering.
    %   U = INITFCM(CLUSTER_N, DATA_N) randomly generates a fuzzy partition
    %   matrix U that is CLUSTER_N by DATA_N, where CLUSTER_N is number of
    %   clusters and DATA_N is number of data points. The summation of each
    %   column of the generated U is equal to unity, as required by fuzzy
    %   c-means clustering.
    %
    %       See also DISTFCM, FCMDEMO, IRISFCM, STEPFCM, FCM.
    
    %   Roger Jang, 12-1-94.
    %   Copyright 1994-2002 The MathWorks, Inc. 
    
    U = rand(cluster_n, data_n);
    col_sum = sum(U);
    U = U./col_sum(ones(cluster_n, 1), :);
end

function [U_new, center, obj_fcn] = stepfcm(data, U, cluster_n, expo)
    %STEPFCM One step in fuzzy c-mean clustering.
    %   [U_NEW, CENTER, ERR] = STEPFCM(DATA, U, CLUSTER_N, EXPO)
    %   performs one iteration of fuzzy c-mean clustering, where
    %
    %   DATA: matrix of data to be clustered. (Each row is a data point.)
    %   U: partition matrix. (U(i,j) is the MF value of data j in cluster j.)
    %   CLUSTER_N: number of clusters.
    %   EXPO: exponent (> 1) for the partition matrix.
    %   U_NEW: new partition matrix.
    %   CENTER: center of clusters. (Each row is a center.)
    %   ERR: objective function for partition U.
    %
    %   Note that the situation of "singularity" (one of the data points is
    %   exactly the same as one of the cluster centers) is not checked.
    %   However, it hardly occurs in practice.
    %
    %       See also DISTFCM, INITFCM, IRISFCM, FCMDEMO, FCM.
    
    %   Copyright 1994-2014 The MathWorks, Inc. 
    
    mf = U.^expo;       % MF matrix after exponential modification
    center = mf*data./(sum(mf,2)*ones(1,size(data,2))); %new center
    dist = distfcm(center, data);       % fill the distance matrix
    obj_fcn = sum(sum((dist.^2).*mf));  % objective function
    tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1
    U_new = tmp./(ones(cluster_n, 1)*sum(tmp));
end