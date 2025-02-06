classdef Mario < PROBLEM
% <2024> <single> <integer/label>
% Play with Mario!
% stageNo ---  1 --- Stage no. (1-3)
% nFrame  --- 60 --- Number of frames per second
% isImage ---  1 --- Show image or not
% isSound ---  1 --- Play sound or not

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        Data;
        map;
        nFrame;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [stageNo,obj.nFrame,isImage,isSound] = obj.ParameterSet(1,60,1,1);
            obj.M = 1;
            if isempty(obj.D); obj.D = 200; end
            obj.lower    = [zeros(1,obj.D/2)+1,zeros(1,obj.D/2)+5];
            obj.upper    = [zeros(1,obj.D/2)+5,zeros(1,obj.D/2)+50];
            obj.encoding = [zeros(1,obj.D/2)+3,zeros(1,obj.D/2)+2];
            obj.map = MAP(18,200);
            load(fullfile(fileparts(mfilename('fullpath')),'images.mat'),'Map');
            obj.Data = Map(stageNo,:);
            if isImage
                obj.map.LoadImage();
            end
            if isSound
                obj.map.LoadSound();
            end
            obj.map.AddTileSprite(obj.Data{1});
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            if nargin > 1
                PopDec(:,end/2+1:end) = cumsum(PopDec(:,end/2+1:end),2);
                PopObj = zeros(size(PopDec,1),1);
                obj.map.DelSprite();
                obj.map.AddTileSprite(obj.Data{2});
                obj.map.AddTileSprite(repmat([12,15,6],size(PopDec,1),1));
                iter = 0;
                while iter<=max(PopDec(:,end)) && obj.map.Action()
                    iter = iter+1;
                    pause(1/obj.nFrame);
                    for i = find(isvalid(obj.map.player))
                        index = find(PopDec(i,end/2+1:end)==iter);
                        if ~isempty(index)
                            switch PopDec(i,index)
                                case 1
                                    obj.map.player(i).moving = 1;
                                case 2
                                    obj.map.player(i).moving = -1;
                                case 3
                                    obj.map.player(i).jumping = 1;
                                case 4
                                    obj.map.player(i).moving = 0;
                                case 5
                                    obj.map.player(i).jumping = 0;
                            end
                        end
                        PopObj(i) = min(PopObj(i),-obj.map.player(i).loc(2));
                    end
                end
            elseif length(obj.map.images{1}) > 1
                obj.map.DelSprite();
                obj.map.AddTileSprite(obj.Data{2});
                obj.map.AddTileSprite([12,15,6]);
                obj.map.Show();
                set(obj.map.imageobj.Parent.Parent,'KeyPressFcn',{@fcn,obj.map,1},'KeyReleaseFcn',{@fcn,obj.map,0});
                while isvalid(obj.map.imageobj.Parent.Parent) && ~obj.map.imageobj.Parent.Parent.BeingDeleted && obj.map.Action()
                    pause(1/obj.nFrame);
                end
            end
        end
    end
end

function fcn(~,event,map,type)
    player = map.player(find(isvalid(map.player),1));
    if ~isempty(player)
        switch event.Key
            case 'uparrow'
                if type
                    player.jumping = 1;
                else
                    player.jumping = 0;
                end
            case 'leftarrow'
                if type
                    player.moving = -1;
                else
                    player.moving = 0;
                end
            case 'rightarrow'
                if type
                    player.moving = 1;
                else
                    player.moving = 0;
                end
        end
    end
end