classdef SPRITE < handle & matlab.mixin.Heterogeneous

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        map;
        images;
        state     = 1;
        loc       = [0 0];
        velocity  = [0 0];
        direction = 1;
        ground    = 1;
        alive     = 1;
    end
    methods
        function obj = SPRITE(loc,images)
            obj.loc    = loc;
            obj.images = images;
        end
        function CollisionBlock(obj)
            obj.loc(2) = obj.loc(2) + obj.velocity(2);
            if any(obj.map.block([ceil(obj.loc(1)/16),ceil((obj.loc(1)+15)/16)],ceil((obj.loc(2)+15)/16)))
                obj.loc(2) = floor((obj.loc(2)-1)/16)*16+1;
            elseif any(obj.map.block([ceil(obj.loc(1)/16),ceil((obj.loc(1)+15)/16)],ceil(obj.loc(2)/16)))
                obj.loc(2) = ceil(obj.loc(2)/16)*16+1;
            end
            obj.loc(1) = obj.loc(1) + obj.velocity(1);
            if any(obj.map.block(ceil(obj.loc(1)/16),[ceil(obj.loc(2)/16),ceil((obj.loc(2)+15)/16)]))
                obj.loc(1)      = ceil(obj.loc(1)/16)*16+1;
                obj.velocity(1) = 0;
            elseif any(obj.map.block(ceil((obj.loc(1)+15)/16),[ceil(obj.loc(2)/16),ceil((obj.loc(2)+15)/16)]))
                obj.loc(1) = floor((obj.loc(1)-1)/16)*16+1;
            end
            if any(obj.map.block(ceil((obj.loc(1)+16)/16),[ceil(obj.loc(2)/16),ceil((obj.loc(2)+15)/16)]))
                obj.velocity(1) = 0;
                obj.ground      = true;
            else
                obj.ground = false;
            end
        end
        function CollisionEnemy(obj)
            if ~isempty(obj.map.enemy) && any(isvalid(obj.map.enemy))
                aliveEnemy = obj.map.enemy(isvalid(obj.map.enemy));
                coll = find(all(abs(obj.loc-cat(1,aliveEnemy.loc))<16,2));
                if ~isempty(coll)
                    if obj.velocity(1) > 0
                        obj.map.Play('kill');
                        [aliveEnemy(coll).alive] = deal(0);
                        obj.velocity(1) = -10;
                    else
                        obj.alive = 0;
                    end
                end
            end
        end
        function CollisionCoin(obj)
            if ~isempty(obj.map.coin) && any(isvalid(obj.map.coin))
                aliveCoin = obj.map.coin(isvalid(obj.map.coin));
                coll = find(all(abs(obj.loc-cat(1,aliveCoin.loc))<16,2));
                if ~isempty(coll)
                    obj.map.Play('coin');
            	    [aliveCoin(coll).alive] = deal(0);
                end
            end
        end
        function DropDie(obj)
            if obj.loc(1) >= size(obj.map.block,1)*16-31
                obj.alive = 0;
            end
        end
    end
end