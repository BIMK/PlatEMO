function Population = FIP_IPG(Problem,PopDec)
    options.lambda = 0.5;              
	options.dim = floor(Problem.N/2);                    
	options.kernel_type = 'rbf';    
	options.gamma = 20;               
	options.T = 10;      
    %寻找边界点和拐点
    SpecialPoints=findSpecialPoints(PopDec);

    Pop_FI=SpecialPoints;
    lastProblem=Problem;
    lastProblem.FE=Problem.FE-Problem.N;
    X=lastProblem.Initialization();
    Y=Problem.Initialization();
    Fx=X.obj;
    Fy=Y.obj;
    [~,~,A]=JDA(X.decs,Fx,Y.decs,Fy,options);


end

function SpcecialPoints=findSpecialPoints(PopDec)
    % 边界点
    [~,index]=min(PopDec);
    boundary_points = PopDec(index,:);
    distance=zeros(length(PopDec),1);
    for i=1:size(PopDec,1)
        distance(i)=point_to_hyperplane(PopDec(i,:),boundary_points);
    end
    [~,index]=max(distance);
    SpcecialPoints=[boundary_points;PopDec(index,:)];
end

function distance = point_to_hyperplane(point, plane_points)
    % point: 要计算距离的点
    % plane_points: 构成超平面的点的坐标
    
    % 使用任意一个点作为超平面方程的参考点
    reference_point = plane_points(1, :);
    
    % 通过其他点计算超平面的法向量
    % 取第一个点作为起点，计算向量
    vectors = plane_points - reference_point;
    
    % 计算法向量
    % normal_vector = null(vectors); % 使用 null 函数求解矩阵的零空间，得到法向量
    % 假设有一个3D空间的点集合

    % 计算点集合的中心
    center = mean(plane_points);
    
    % 将点集合减去中心，得到零均值的数据
    centered_points = plane_points - center;
    
    % 计算协方差矩阵
    covariance_matrix = cov(centered_points);
    
    % 对协方差矩阵进行特征值分解
    [eigenvectors, eigenvalues] = eig(covariance_matrix);
    
    % 法向量就是最小特征值对应的特征向量
    [~, min_index] = min(diag(eigenvalues));
    normal_vector = eigenvectors(:, min_index);


    % 超平面方程的常数项
    % 根据 Ax + By + Cz + D = 0 的形式，D = -Ax - By - Cz
    D = -normal_vector' * reference_point';
    
    % 计算点到超平面的距离
    distance = abs(normal_vector' * point' + D) / norm(normal_vector);
end


