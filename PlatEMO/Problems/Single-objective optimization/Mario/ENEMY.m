classdef ENEMY < SPRITE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function obj = ENEMY(loc,images)
            obj = obj@SPRITE(loc,images);
            obj.velocity = [0 -1];
        end
        function Action(obj)
            if obj.alive > 0
                if any(obj.map.block([ceil(obj.loc(1)/16),ceil((obj.loc(1)+15)/16)],ceil((obj.loc(2)+16)/16)))
                    obj.velocity(2) = -abs(obj.velocity(2));
                elseif any(obj.map.block([ceil(obj.loc(1)/16),ceil((obj.loc(1)+15)/16)],ceil((obj.loc(2)-1)/16)))
                    obj.velocity(2) = abs(obj.velocity(2));
                end
                if ~obj.ground
                    obj.velocity(1) = min(16,obj.velocity(1)+1);
                end
                obj.CollisionBlock();
                obj.DropDie();
            elseif obj.alive > -10
                obj.alive = obj.alive-1;
            else
                delete(obj);
                return;
            end
            if obj.alive <= 0
                obj.state = 3;
            else
                obj.state = mod(obj.state,2)+0.2;
            end
        end
    end
end