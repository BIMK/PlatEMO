function k = MRP(obj_g,obj_g_1)
% 划分子区域
    u=max(obj_g);
    l=min(obj_g);
    delta_f=(obj_g-obj_g_1)./(u-l);
    mu=mean(abs(delta_f));
    delta=mean(mean(delta_f-repmat(mu,size(delta_f,1),1)));
    M=size(obj_g,2);
    k1=M+1;
    k2=M*2;
    k=ceil(k1+delta*(k2-k1));
    fprintf("k=%d\n",k);
    [W,~]=UniformPoint(k,size(obj_g,2),"MUD");
% 个体关联
    distance=zeros(size(obj_g,1),k);
    for i = 1:size(distance,1)
        for j = 1:size(distance,2)
            distance(i,j)=Point2VectorDistance(obj_g(i,:),W(j,:));
        end
    end
    
    [~,h]=min(distance,[],2);
    % 先验证分区是否正确
    % figure
    % gscatter(obj_g(:,1),obj_g(:,2),h);
    % hold on
    % for i = 1:size(W, 1)
    %     direction_vector = W(i,:);
    %     plot([0, direction_vector(1)], [0, direction_vector(2)], 'LineWidth', 2);
    % end

% 确定进化步长
    

end

function distance = Point2VectorDistance(point, vector)
    % point: 点的坐标 [x, y, z]
    % vector: 空间向量的终点坐标 [x, y, z]
    
    % 计算点到向量的投影点
    projection_point = (dot(point, vector) / dot(vector, vector)) * vector;
    
    % 计算点到投影点的距离
    distance = norm(point - projection_point);
end