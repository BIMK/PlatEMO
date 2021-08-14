function net = GNGUpdate(oSignals,net)
% Update GNG

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
    
    %% Parameters 
    [nSig,D]=size(oSignals);
    maxAge = net.maxAge; 
    lambda = net.lambda;   
    epsilon_a = net.epsilon_a;    
    epsilon_nb = net.epsilon_nb;   
    alpha = net.alpha;           
    delta = net.delta;           
    maxHP = net.maxHP;
    maxIter = net.maxIter;
    Node = net.Node;
    Err  = net.Err;   
    edge = net.edge;  
    age  = net.age;      
    hp  = net.hp;

    
    %% Randamize Input Signal
    ran = randperm(nSig);
    oSignals = oSignals(ran,:);   
   
    %% Input Signal Normalization [0-1]
    oMax = max(oSignals,[],1);
    oMin = min(oSignals,[],1);
    oRange = oMax-oMin;
    Signals = (oSignals - repmat(oMin,nSig,1))./repmat(oRange,nSig,1);

    %% Update GNG   
    for nitr = 1:maxIter      
        % Step 0: Initialization. Start with two neural units (nodes) selected from input data
        if size(Node,1) <= 2
            Ni = 2;
            Xmin = min(Signals,[],1);
            Xmax = max(Signals,[],1);
            for i = 1:Ni
                Node(i,:) = unifrnd(Xmin, Xmax);
            end     
            Err = [0; 0];
            edge = zeros(2,2);  
            age  = zeros(2,2);          
            hp = ones(1,2).*maxHP;
        end
                       
        for numSig = 1:nSig

            % Step 1: Input one signal
            pattern = Signals(numSig,:);

            % Step 2: Find the two nearest nodes ra and rb to new signal
            d = pdist2(pattern, Node);
            [~, SortOrder] = sort(d);
            ra = SortOrder(1);
            rb = SortOrder(2);

            % Step 2.5: change HP
            hp = hp-1;
            hp(ra) = maxHP;
            hp(rb) = hp(rb)+1;

            % Steps 3: Increment the age of all edges emanating from ra
            age(ra, :) = age(ra, :) + 1;
            age(:, ra) = age(:, ra) + 1;

            % Step 4: Add the squared distance to a local error counter variable
            Err(ra) = Err(ra) + d(ra)^2;    

            % Step 5: Move ra and its topological neighbors towards singal         
            Node(ra,:) = Node(ra,:) + epsilon_a * (pattern - Node(ra,:));
            for j = find(edge(ra,:)==1)
                Node(j,:) = Node(j,:) + epsilon_nb * (pattern - Node(j,:));% for Neighbor nodes which are connecting to ra.      
            end

            % Step 6:
            % If ra and rb are connected by an edge, set the age of this edge to zero.
            % If such an edge does not exist, create it.
            edge(ra,rb) = 1;
            edge(rb,ra) = 1;
            age(ra,rb) = 0;
            age(rb,ra) = 0;

            % Step 7(1):
            % Remove edges from node if age>maxAge.
            edge(age>maxAge) = 0;

            % Step 7(2):
            % Remove dead node and their edges.       
            DeadNodes = (hp<=0);
            edge(DeadNodes, :) = [];
            edge(:, DeadNodes) = [];
            age(DeadNodes, :) = [];
            age(:, DeadNodes) = [];
            Node(DeadNodes, :) = [];
            Err(DeadNodes) = [];
            hp(DeadNodes) = [];

            % Step 8: Node Insertion Procedure.
            % rnew: new node
            % r1max: node which has maximum accumulated error
            % r2max: neighbor node of r1max
            if mod(numSig, lambda) == 0 && net.maxNode > size(Node,1)
                [~, r1max] = max(Err);
                [~, r2max] = max(edge(:,r1max).*Err);
                rnew = size(Node,1) + 1;
                Node(rnew,:) = (Node(r1max,:) + Node(r2max,:))/2;   
                edge(r1max,r2max) = 0;  
                edge(r2max,r1max) = 0;
                edge(r1max,rnew) = 1;  
                edge(rnew,r1max) = 1;
                edge(rnew,r2max) = 1;  
                edge(r2max,rnew) = 1;
                age(rnew,:) = 0;   
                age(:,rnew) = 0;
                Err(r1max) = alpha * Err(r1max); 
                Err(r2max) = alpha * Err(r2max);
                Err(rnew) = Err(r1max);  
                hp(rnew) = maxHP;
            end

            % Step 9: Decrease the error of all units.
            Err = delta*Err;

        end     
    end     
    
    %% Output net
    net.Node = Node;
    net.Err = Err;
    net.edge = edge;
    net.age = age;
    net.hp = hp;
    
   
    %% Expansion (Algorithm 3)
    % Accociate signal to its closest node
    Distance = pdist2(Node,Signals);
    [~,pi] = min(Distance,[],1); 
    % Node Label based on edge (detect sub-network)
    connection = graph(edge ~= 0);
    NetLabel = conncomp(connection);
    % Data Label based on Node 
    DataLable = NetLabel(pi);
    % Expand each sub-network 
    for i=1:max(NetLabel)
        subData = find(DataLable == i);
        subNet = find(NetLabel == i);       
        if length(subData)<=1 || length(subNet)<=1
            continue;
        end        
        DataMax = max(Signals(subData,:),[],1);
        DataMin = min(Signals(subData,:),[],1);
        NetMax = max(Node(subNet,:),[],1);
        NetMin = min(Node(subNet,:),[],1);
        DataRange = DataMax-DataMin;
        NetRange = NetMax-NetMin;
        Ratio = DataRange./NetRange;
        for k = 1:D
            for j = subNet
                Node(j,k)= (Node(j,k)- NetMin(k)).*Ratio(k)+DataMin(k);
            end
        end
    end
    net.NodeS = Node;   
end