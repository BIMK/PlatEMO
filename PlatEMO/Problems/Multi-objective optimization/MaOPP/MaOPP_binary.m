classdef MaOPP_binary < PROBLEM
% <many> <binary> <large/none> <expensive/none>
% Many-objective pathfinding problem based on binary encoding
% xmax                    ---  10 --- xmax
% ymax                    ---  10 --- ymax
% obstacleValue           ---   0 --- obstacleValue
% nh                      ---   1 --- nh
% neighbourhood           ---   2 --- neighbourhood
% backtracking            ---   0 --- backtracking
% allowObstaclesOnPath    ---   1 --- allowObstaclesOnPath
% overheadVariablesFactor --- 1.5 --- overheadVariablesFactor

%------------------------------- Reference --------------------------------
% J. Weise and S. Mostaghim, A scalable many-objective pathfinding
% benchmark suite, IEEE Transactions on Evolutionary Computation, 2022,
% 26(1): 188-194.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The most recent reference fronts and sets can be found at 
% https://ci.ovgu.de/Publications/Scalable+Many_Objective+Pathfinding+Benchmark+Suite-p-910.html
% In case of any question, please contact me at 
% jens.weise@ovgu.de

    properties(Access = protected)
        x_max;          % size [1,x_max] in x direction
        y_max;          % size [1,y_max] in y direction
        delay;          % expected delay
        nh;             % elevation function {1,2,3,M} -> 1 = 1, 2 = 2, 3 = 3, 4 = M
        neighborhood;   % 2^k, k {2,3}
        backtracking;   % True, False
        obstacle;       %{0,1,2} -> 0 = No, 1 = CH, 2 = LA
        vmax_high;
        vmax_medium;
        vmax_low;
        v_max;
        elevation;
        allowObstaclesOnPath;
        upperBoundsForObjectives;
        overheadVariablesFactor;
        numerOfBits;
    end

    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.x_max, obj.y_max, obstacleValue, obj.nh, obj.neighborhood, obj.backtracking, obj.allowObstaclesOnPath, obj.overheadVariablesFactor] = obj.ParameterSet(10,10,0,1,2,0,1,1);
            if isempty(obj.vmax_high); obj.vmax_high = 130; end
            if isempty(obj.vmax_medium); obj.vmax_medium = 100; end
            if isempty(obj.vmax_low); obj.vmax_low = 50; end
            obj.M = 5;
            obj.setBinaryEncoding();
            obj.v_max     = zeros(obj.x_max, obj.y_max);  % matrix of shape (x_max,y_max), matrix with different v values for each cell
            obj.elevation = zeros(obj.x_max, obj.y_max);  % matrix of shape (x_max,y_max)
            for x = 1 : obj.x_max
                for y = 1 : obj.y_max
                    obj.obstacle(x,y) = obstacleValue;
                    w_x_y = max(sin(x-1), cos(y-1));
                    if w_x_y > 0.9
                        obj.v_max(x, y) = obj.vmax_high;
                    elseif w_x_y < -0.4
                        obj.v_max(x, y) = obj.vmax_low;
                    else
                        obj.v_max(x, y) = obj.vmax_medium;
                    end
                    % velocity matrix obstacles for CH
                    if obj.obstacle(x,y) == 1 && (sign(sin(pi/2 + pi*x)) + sign(sin(pi/2 + pi*y)) - 2* (heaviside(x-obj.x_max+0.5)- heaviside(x-obj.x_max-0.5))* (heaviside(y-obj.y_max+0.5)- heaviside(y-obj.y_max-0.5)) == 2)
                        obj.v_max(x, y) = 0;
                    end
                    % velocity matrix obstacles for LA
                    if obj.obstacle(x,y) == 2 && (x-1-obj.x_max/2)^2 + (y-1-obj.y_max/2)^2 - (0.25*obj.x_max)^2 <0 %Q: using radius of x_max/4 from paper page 4 line 22 A: Yes.
                        obj.v_max(x, y) = 0;
                    end
                    xs = map(x,1,obj.x_max+1,-3,3);
                    ys = map(y,1,obj.y_max+1,-3,3);
                    if obj.nh == 1
                        obj.elevation(x,y) = 5*exp(-(xs-1.5)^2-(ys+1.5)^2);
                    elseif obj.nh == 2
                        obj.elevation(x,y) = 5*exp(-(xs+1.5)^2-(ys+1.5)^2) + 5*exp(-(xs-1.5)^2-(ys-1.5)^2);
                    elseif obj.nh == 3
                        obj.elevation(x,y) = 5*exp(-(xs+1.5)^2-(ys+1.5)^2) + 5*exp(-(xs-1.5)^2-(ys-1.5)^2) + 5*exp(-(xs-1.5)^2-(ys+1.5)^2);
                    elseif obj.nh == 4
                        obj.elevation(x,y) = 3*(1-xs)^2*exp(-(xs^2)-(ys+1)^2)-10*exp(-xs^2-ys^2)*(-xs^3+xs/5-ys^5)-1/3*exp(-(xs+1)^2-ys^2);
                    end
                end
            end
            manhattanSteps = (obj.x_max + obj.y_max - 2);
            obj.upperBoundsForObjectives = [manhattanSteps*1.5, manhattanSteps*1.5, 1.5*5*obj.nh, 1.5*manhattanSteps/50, 0.5*manhattanSteps*pi/2];
        end
        function setBinaryEncoding(obj)
            if obj.neighborhood == 2 && obj.backtracking == 0
                obj.D = 1*(obj.x_max + obj.y_max - 2)*obj.overheadVariablesFactor;
                obj.numerOfBits = 1;
            elseif obj.neighborhood == 2 && obj.backtracking == 1
                obj.D = 2*(obj.x_max + obj.y_max - 2)*obj.overheadVariablesFactor;
                obj.numerOfBits = 2;
            elseif obj.neighborhood == 3 && obj.backtracking == 0
                obj.D = 2*(obj.x_max + obj.y_max - 2)*obj.overheadVariablesFactor;
                obj.numerOfBits = 2;
            elseif obj.neighborhood == 3 && obj.backtracking == 1
                obj.D = 3*(obj.x_max + obj.y_max - 2)*obj.overheadVariablesFactor;
                obj.numerOfBits = 3;
            end
            obj.encoding = 4 + zeros(1,obj.D);
        end
        function [x_coords, y_coords,D] = decodePath(obj,PopDec,x_max,y_max)
            D = size(PopDec,2);
            n = obj.numerOfBits;
            for i = 0 : (D/n)-1
               temp = PopDec(:,[n*i+1:n*i+n]);
               realPopDec(:,i+1) = b2d(temp); 
            end
            realPopDec = 1.0*realPopDec/(2^n);
            [x_coords, y_coords, D] = decodePathWithRealValues(obj,realPopDec,x_max,y_max);
            D = size(realPopDec,2);
        end
        function [x_coords, y_coords, D] = decodePathWithRealValues(obj,PopDec,x_max,y_max)
            [N,D] = size(PopDec);
            % convert the array into a sequence of x,y coordinates
            x_coords = zeros(N,D+1);
            y_coords = zeros(N,D+1);
            x_coords(:,1) = 1;
            y_coords(:,1) = 1;

            for i = 1 : N
                x_current = 1;
                y_current = 1;
                for d = 1 : D
                    possibleNeighours      = obj.getPossNeighbours(x_current,y_current,obj.backtracking,obj.neighborhood,obj.allowObstaclesOnPath);
                    noOfPossibleNeighbours = size(possibleNeighours,1);
                    neighbourToTake        = noOfPossibleNeighbours;
                    for a = 1 : noOfPossibleNeighbours
                        if PopDec(i,d) >= (a-1)/noOfPossibleNeighbours && PopDec(i,d) < (a)/noOfPossibleNeighbours
                            neighbourToTake = a;
                            break;
                        end
                    end
                    if neighbourToTake == 0 % Path ends if there are no neighbours available
                        break;
                    end
                    addvals   = possibleNeighours(neighbourToTake,:);
                    xAdd      = addvals(1);
                    yAdd      = addvals(2);
                    x_current = x_current + xAdd;
                    y_current = y_current + yAdd;
                    x_coords(i,d+1) = x_current;
                    y_coords(i,d+1) = y_current;
                    if x_current == obj.x_max && y_current == obj.y_max
                        break;
                    end
                end
            end
        end
        function neighbours = getPossNeighbours(obj,xCurr,yCurr,backtracking,neighborhood,allowObstaclesOnPath)
            e  = [1 0];
            se = [1 1];
            s  = [0 1];
            sw = [-1 1];
            w  = [-1 0];
            nw = [-1 -1];
            n  = [0 -1];
            ne = [1 -1];
            currentCoordinates = [xCurr yCurr];
            if backtracking == 0
                if neighborhood == 2
                    neighbours = [e;s];
                else
                    neighbours=[e;se;s];
                end
            else
                if neighborhood == 2
                    neighbours = [e;s;w;n];
                else
                    neighbours = [e;se;s;sw;w;nw;n;ne];
                end
            end
            % should we remove the outliers?
            % it is closer to the original implementation
            nextCoordinates = neighbours+currentCoordinates;
            out = nextCoordinates(:,1)<1 | nextCoordinates(:,1)>obj.x_max | nextCoordinates(:,2)<1 | nextCoordinates(:,2)>obj.y_max;
            neighbours(out,:) = [];
            if allowObstaclesOnPath == 0
                nextCoordinates = neighbours+currentCoordinates;
                obstancles = getVmaxValues(nextCoordinates(:,2),nextCoordinates(:,1),obj.v_max) == 0;
                neighbours(obstancles,:) = [];
            end
        end
        %% separate the calculations here to be able to compute the feasibility rate as well (with a lot of extra runtime, but whatever..)
        function [objectives,distToTarget,fObstacle] = objectiveValuesPath(obj,PopDec)
            % read out parameters
            N = size(PopDec,1);
            [x_coords, y_coords, D] = obj.decodePath(PopDec,obj.x_max,obj.y_max);
            % for entire population (one column) iterate over nodes (columns)
            f1 = zeros(N,1); % euclidean distance
            f2 = zeros(N,1); % expected delays
            f3 = zeros(N,1); % elevations
            f4 = zeros(N,1); % traveling time
            f5 = zeros(N,1); % smoothness
            fObstacle = zeros(N,1);
            distances = zeros(N,D); % stores distances that were calculated for f1 for later use in f5
            c = [];
            for d = 1 : D
                x1        = x_coords(:,d);
                x2        = x_coords(:,d+1);
                y1        = y_coords(:,d);
                y2        = y_coords(:,d+1);
                hasEnded  = (x2~=0&x1~=0);
                x2(x2==0) = 1;
                y2(y2==0) = 1;
                x1(x1==0) = 1;
                y1(y1==0) = 1;
                
                vMaxValues1 = getVmaxValues(x1,y1,obj.v_max);
                vMaxValues2 = getVmaxValues(x2,y2,obj.v_max);
                
                if obj.allowObstaclesOnPath==1
                    obstancles = vMaxValues1==0|vMaxValues2==0;
                    fObstacle  = fObstacle+obstancles;
                else
                    obstancles = zeros(N,1);
                    fObstacle  = fObstacle+obstancles;
                end

                %---Objective 1: euclidean distance
                distances(:,d) = sqrt((x2-x1).^2+(y2-y1).^2);
                res = (distances(:,d)).*hasEnded;
                res(res>sqrt(2)) = 0;
                f1 = f1 + res;

                %---Objective 2: delays
                inds      = vMaxValues1 ~= vMaxValues2;
                del       = inds*2;
                inds      = vMaxValues1 == vMaxValues2;
                loc       = find(inds==1);
                inds      = zeros(N,1);
                inds(loc) = vMaxValues1(loc);
                inds50    = inds == 50;
                del       = del + inds50*3;
                inds100   = inds == 100;
                del       = del + inds100*1;
                inds0     = del == 0;
                del       = del + inds0*0.2;
                res       = del.*hasEnded;
                f2        = f2 + res;

                %---Objective 3: elevations
                h1   = getVmaxValues(x1,y1,obj.elevation);
                h2   = getVmaxValues(x2,y2,obj.elevation);
                inds = h1 < h2;
                res  = (inds.*(h2-h1)).*hasEnded;
                f3   = f3 + res;

                %---Objective 4: traveling time
                res = (2*distances(:,d)./(vMaxValues1+vMaxValues2));
                res(isinf(res)|isnan(res)) = 0;
                if any(isnan(res))
                    errr;
                end
                res = res.*hasEnded;
                f4  = f4 + res;

                %---Objective 5: smoothness
                if d>1
                    x0       = x_coords(:,d-1);
                    y0       = y_coords(:,d-1);
                    x0(x0<1) = 1;
                    y0(y0<1) = 1;
                    
                    test = transpose(dot([(x0-x1).'; (y0-y1).'],[(x1-x2).'; (y1-y2).'],1))./(distances(:,d-1).*distances(:,d));
                    test(isnan(test)) = 0;
                    try
                    	assert(all(test <= 1) && all(test >= -1))
                    catch
                    end
                    res = (acos(transpose(dot([(x0-x1).'; (y0-y1).'],[(x1-x2).'; (y1-y2).'],1))./(distances(:,d-1).*distances(:,d)))).*hasEnded;
                    res(isnan(res) )= 0;
                    f5 = f5 + res;
                end
            end
            objectives = [f1,f2,f3,f4,f5];
            try
                last = [];
                for i = 1 : N
                    last = [last find(x_coords(i,:)~=0,1,'last')];
                end
                last = last';
                lxc  = [];
                for k = 1 : N
                    lxc(k,1) = x_coords(k,last(k));
                    lxc(k,2) = y_coords(k,last(k));
                end
            catch
            end
            target       = [repmat(obj.x_max,N,1) repmat(obj.y_max,N,1)];
            distToTarget = abs(target(:,1)-lxc(:,1))+abs(target(:,2)-lxc(:,2));
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            [objectives,distToTarget,fObstacle] = objectiveValuesPath(obj,varargin{1});
            objectives(:,1)                     = objectives(:,1) + 2*distToTarget + fObstacle;
            Population = SOLUTION(varargin{1},objectives,distToTarget+fObstacle,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        % Here, the pre-computed fronts are read. If it is not available,
        % an educated-guess Nadir point is computed. Please check the
        % website above in the header to see, if there are new fronts
        % available.
        function R = GetOptimum(obj,N)
            obs = '';
            if obj.obstacle == 0
                obs = 'NO';
            elseif obj.obstacle == 1
                obs = 'CH';
            elseif obj.obstacle == 2
                obs = 'LA';
            end
            sizeX = obj.x_max;
            sizeY = obj.y_max;
            ele   = '';
            if obj.nh == 1
                ele = '1';
            elseif obj.nh == 2
                ele = '2';
            elseif obj.nh == 3
                ele = '3';
            elseif obj.nh == 4
                ele = 'M';
            end
            k = obj.neighborhood;
            if obj.backtracking
                bt = 'T';
            else
                bt = 'F';
            end
            problem = sprintf('ASLETISMAC_%s_X%i_Y%i_P%s_K%i_B%s',obs,sizeX,sizeY,ele,k,bt);
            try
                CallStack = dbstack('-completenames');
                load(fullfile(fileparts(CallStack(1).file),'PF.mat'),problem);
                order = [2 5 1 3 4];
                eval(['R=',problem,';']);
                R      = R(:,order);
                R(:,5) = deg2rad(R(:,5));
            catch
                R = repmat(obj.upperBoundsForObjectives,N,1);
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end

function out= map(x,in_min,in_max,out_min,out_max)
    out=(x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
end

function vList = getVmaxValues(x1,y1,aMatrix)
    len = size(x1,1);
    vList = zeros(len,1);
    for i = 1:len
        vList(i) = aMatrix(x1(i),y1(i));
    end
end

% Converts binary to decimal (MSB 0)!
function dec = b2d(vector)
    l   = size(vector);
    dec = ones(l(1),1);
    for r = 1 : l(1)
        dTemp = 0;
        for i = 1 : l(2)
            dTemp = dTemp + vector(r,i)*2^(i-1);
        end
        dec(r,1) = dTemp;
    end
end