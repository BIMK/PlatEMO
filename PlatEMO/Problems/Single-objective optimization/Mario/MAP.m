classdef MAP < handle

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        images = {};
        sounds = struct();
        background;
        imageobj;
        block;
        player;
        enemy;
        coin;
    end
    methods
        function obj = MAP(xsize,ysize)
            obj.background = repmat(reshape(uint8([115,146,241]),1,1,3),xsize*16,ysize*16,1);
            obj.block      = false(xsize,ysize);
            obj.block([1,end],:) = true;
            obj.block(:,[1,end]) = true;
            [obj.images{1:21}]   = deal(0);
            obj.images{8}        = zeros(32,32,3);
            obj.images{9}        = zeros(16,32,3);
        end
        function AddTileSprite(obj,data)
            pixels = (data(:,2:3)-1)*16+1;
            for i = 1 : size(data,1)
                if data(i,1) < 10
                    tile = obj.images{data(i,1)};
                    obj.background(pixels(i,1):pixels(i,1)+size(tile,1)-1,pixels(i,2):pixels(i,2)+size(tile,2)-1,:) = tile;
                    obj.block(data(i,2):data(i,2)+ceil(size(tile,1)/16)-1,data(i,3):data(i,3)+ceil(size(tile,2)/16)-1) = data(i,1)>6;
                elseif data(i,1) == 10
                    obj.coin = [obj.coin,COIN(pixels(i,1:2),obj.images([10:12,11]))];
                    obj.coin(end).map = obj;
                elseif data(i,1) == 11
                    obj.enemy = [obj.enemy,ENEMY(pixels(i,1:2),obj.images(13:15))];
                    obj.enemy(end).map = obj;
                elseif data(i,1) == 12
                    obj.player = [obj.player,PLAYER(pixels(i,1:2),obj.images([16:18,17,19:21]))];
                    obj.player(end).map = obj;
                end
            end
        end
        function DelSprite(obj)
            delete([obj.player,obj.enemy,obj.coin]);
            obj.player = [];
            obj.enemy  = [];
            obj.coin   = [];
            if isfield(obj.sounds,'bgm')
                obj.sounds.bgm{1}.StopFcn = @(h,~)play(h);
                play(obj.sounds.bgm{1});
            end
        end
        function running = Action(obj)
            allSprite = [obj.player,obj.enemy,obj.coin];
            for sprite = allSprite(isvalid(allSprite))
            	sprite.Action();
            end
            running = any(isvalid(obj.player));
            if running
                obj.Show();
            elseif isfield(obj.sounds,'bgm')
                obj.sounds.bgm{1}.StopFcn = [];
                stop(obj.sounds.bgm{1});
            end
        end
        function LoadImage(obj)
            if length(obj.images{1}) < 2
                load(fullfile(fileparts(mfilename('fullpath')),'images.mat'),'images');
                obj.images = images;
            end
        end
        function Show(obj)
            if length(obj.images{1}) > 1
                img = obj.background;
                allSprite = [obj.player,obj.enemy,obj.coin];
                for sprite = allSprite(isvalid(allSprite))
                    if sprite.direction > 0
                        spriteimg = sprite.images{max(1,round(sprite.state))};
                    else
                        spriteimg = fliplr(sprite.images{max(1,round(sprite.state))});
                    end
                    loc   = {sprite.loc(1):sprite.loc(1)+size(spriteimg,1)-1,sprite.loc(2):sprite.loc(2)+size(spriteimg,2)-1};
                    im    = img(loc{1},loc{2},:);
                    index = repmat(any(spriteimg<185,3),1,1,3);
                    im(index) = spriteimg(index);
                    img(loc{1},loc{2},:) = im;
                end
                alivePlayer  = obj.player(isvalid(obj.player));
                locs         = cat(1,alivePlayer.loc);
                [~,leader]   = max(locs(:,2));
                cameraLoc(1) = min(size(obj.background,1)-271,max(17,alivePlayer(leader).loc(1)-128));
                cameraLoc(2) = min(size(obj.background,2)-271,max(17,alivePlayer(leader).loc(2)-120));
                if isempty(obj.imageobj) || ~isvalid(obj.imageobj)
                    fig = figure('MenuBar','none','ToolBar','none','NumberTitle','off','Position',[500 100 500 500],'Name','Running ...');
                    movegui(fig,'center');
                    obj.imageobj = imshow(img(cameraLoc(1):cameraLoc(1)+255,cameraLoc(2):cameraLoc(2)+255,:),'Parent',axes(fig,'Position',[0 0 1 1]));
                else
                    obj.imageobj.CData = img(cameraLoc(1):cameraLoc(1)+255,cameraLoc(2):cameraLoc(2)+255,:);
                end
                set(obj.imageobj.Parent.Parent,'Name',sprintf('Distance: %d',alivePlayer(leader).loc(2)));
            end
        end
        function LoadSound(obj)
            if ~isfield(obj.sounds,'bgm')
                Cut = @(y)y(find(any(abs(y)>1e-1,2),1):find(any(abs(y)>1e-1,2),1,'last'),:);
                Obj = @(y,fs){audioplayer(y,fs),y,fs};
                Y   = audioread(fullfile(fileparts(mfilename('fullpath')),'sounds.ogg'));
                obj.sounds.bgm  = Obj(Cut(Y(1:842921)),22050);
                obj.sounds.coin = Obj(Cut(Y(842922:868105)),44100);
                obj.sounds.jump = Obj(Cut(Y(868106:873384)),11025);
                obj.sounds.kill = Obj(Cut(Y(873385:878464)),22050);
                obj.sounds.over = Obj(Cut(Y(878465:920947)),22050);
            end
        end
        function Play(obj,str)
            if isfield(obj.sounds,str)
                if isplaying(obj.sounds.(str){1})
                    sound(obj.sounds.(str){2},obj.sounds.(str){3});
                else
                    play(obj.sounds.(str){1});
                end
            end
        end
    end
end