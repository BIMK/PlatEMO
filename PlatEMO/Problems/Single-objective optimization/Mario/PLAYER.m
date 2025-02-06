classdef PLAYER < SPRITE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        moving  = 0;
        jumping = 0;
    end
    methods
        function Action(obj)
            if obj.alive > 0
                if obj.moving > 0
                    obj.velocity(2) = min(5,obj.velocity(2)+1);
                    obj.direction   = 1;
                elseif obj.moving < 0
                    obj.velocity(2) = max(-5,obj.velocity(2)-1);
                    obj.direction   = -1;
                elseif obj.velocity(2) > 0
                    obj.velocity(2) = max(0,obj.velocity(2)-1);
                    obj.direction   = 1;
                elseif obj.velocity(2) < 0
                    obj.velocity(2) = min(0,obj.velocity(2)+1);
                    obj.direction   = -1;
                end
                if obj.jumping && obj.ground
                    obj.velocity(1) = -12;
                    obj.map.Play('jump');
                elseif ~obj.ground
                    obj.velocity(1) = min(16,obj.velocity(1)+1);
                end
                obj.CollisionBlock();
                obj.DropDie();
                obj.CollisionCoin();
                obj.CollisionEnemy()
            elseif obj.alive == 0
                obj.alive    = obj.alive-1;
                obj.velocity = [-10 0];
                obj.map.Play('over');
            elseif obj.alive > -20
                obj.alive       = obj.alive-1;
                obj.loc(1)      = obj.loc(1) + obj.velocity(1);
                obj.velocity(1) = obj.velocity(1)+1;
            else
                delete(obj);
                return;
            end
            if obj.alive <= 0
                obj.state = 7;
            elseif ~obj.ground
                obj.state = 6;
            elseif obj.velocity(2) ~= 0
                obj.state = mod(obj.state-1,4)+1.3;
            else
                obj.state = 1;
            end
        end
    end
end