classdef EnvRL < handle
% Determine environemental selection via reinforcement learning

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        name;   % environemtal selection name
        problem;
        NumEnv; % candidate environmental selection number
        NumState; % state number
        Envs;
        ppo;

        % Basic parameter
        ArcSize;
        Zmin;
        Zmax;
        Cmax;
        Zmin1;
        Zmax1;
        Cmax1;
        Center;
        CenterI;
        CenterII;
        %CenterIII; %empty set
        CenterIV;

        NRI;
        NRII;
        NRIV;
    end

    methods
        function obj = EnvRL(problem, Archive, MaxGen)
            obj.name    = 'EnvRL';
            obj.problem = problem;
            obj.NumEnv  = 3;
            obj.NumState= 10;
            obj.Envs    = {EnvCDP(), EnvEpsilon(), EnvMOP()};
            obj.ppo     = PPO(obj.NumEnv, obj.NumState, MaxGen);
            obj.ArcSize = length(Archive);
            obj.Zmin    = min(Archive.objs, [], 1);
            obj.Zmax    = max(Archive.objs, [], 1);
            obj.Cmax    = max(max(Archive.cons, 0), [], 1);
            obj.Cmax    = max(obj.Cmax, 1e-6);

            % Normalization
            SConsN1 = max(Archive.cons, 0) ./ repmat(obj.Cmax, obj.ArcSize, 1);
            SConsN1 = sum(SConsN1, 2);
            SObjsN1 = (Archive.objs-repmat(obj.Zmin, obj.ArcSize, 1)) ./ repmat(obj.Zmax-obj.Zmin, obj.ArcSize, 1);
            SObjsN1 = sum(SObjsN1, 2);

            obj.Zmin1 = min(SObjsN1);
            obj.Zmax1 = max(SObjsN1);
            obj.Cmax1 = max(SConsN1);

            % Gap solutinon
            CValue     = min(SConsN1);
            OValue     = min(SObjsN1(SConsN1==CValue));
            obj.Center = [OValue, CValue];

            % Re-normalization
            [Sobjs, Scons] = obj.Norm(SObjsN1,SConsN1);

            %% Centroid point
            % I Quadrant: +obj, +con
            IndexI = Sobjs >= 0 & Scons > 0;
            if any(IndexI)
                obj.CenterI = mean([Sobjs(IndexI), Scons(IndexI)],1);
            else
                obj.CenterI = [1,1];
            end

            % II Quadrant: -obj, +con
            IndexII = Sobjs < 0 & Scons > 0;
            if any(IndexII)
                obj.CenterII = mean([Sobjs(IndexII), Scons(IndexII)],1);
            else
                obj.CenterII = [0,0];
            end

            % III Quadrant: -obj, -con
            % empty set

            % IV Quadrant: +obj, -con
            IndexIV = Sobjs > 0 & Scons <= 0;
            if any(IndexIV)
                obj.CenterIV = min([Sobjs(IndexIV), Scons(IndexIV)],[],1);
            else
                obj.CenterIV = [1, 0];
            end
        end

        function [obj, Population, action] = do(obj, Population, VAR, Archive)
      	    % Choose environmental selection
            [state, reward]   = obj.GetStateReward(Archive);            
            [obj.ppo, action] = obj.ppo.GetAction(state, reward);
            % Population, N
            if action == 1
                Population = obj.Envs{1}.do(Population, obj.problem.N);
            elseif action == 2
                Population = obj.Envs{2}.do(Population, obj.problem.N, VAR);
            else
                Population = obj.Envs{3}.do(Population, obj.problem.N, VAR);
            end
        end


        function [state, reward] = GetStateReward(obj, Archive)
            %% Calculate Reward
            % Normalization
            SConsN1 = max(Archive.cons, 0) ./ repmat(obj.Cmax, obj.ArcSize, 1);
            SConsN1 = sum(SConsN1, 2);
            SObjsN1 = (Archive.objs-repmat(obj.Zmin, obj.ArcSize, 1)) ./ repmat(obj.Zmax-obj.Zmin, obj.ArcSize, 1);
            SObjsN1 = sum(SObjsN1, 2);

            [Sobjs, Scons] = obj.Norm(SObjsN1, SConsN1);

            % I Quadrant: +obj, +con
            IndexI = Sobjs >= 0 & Scons > 0;
            if any(IndexI)
                CenterI = mean([Sobjs(IndexI), Scons(IndexI)],1);
            else
                CenterI = [1, 1];
            end
            Reward1 = (obj.CenterI(2) - CenterI(2)) + (obj.CenterI(1) - CenterI(1));
            
            % II Quadrant: -obj, +con
            IndexII = Sobjs < 0 & Scons > 0;
            if any(IndexII)
                CenterII = mean([Sobjs(IndexII), Scons(IndexII)],1);
            else
                CenterII = [0,0];
            end
            Reward2 = (obj.CenterII(2) - CenterII(2));
            
            % III Quadrant: -obj, -con
            IndexIII = Sobjs <= 0 & Scons <= 0;
            IndexIII(Sobjs==0 & Scons==0) = false;
            if any(IndexIII)
                CenterIII = min([Sobjs(IndexIII), Scons(IndexIII)],[],1);
            else
                CenterIII = [0,0];
            end
            Reward3 = -(CenterIII(1))-(CenterIII(2));

            % IV Quadrant: +obj, -con
            IndexIV = Sobjs > 0 & Scons <= 0;
            if any(IndexIV)
                CenterIV = min([Sobjs(IndexIV), Scons(IndexIV)],[],1);
            else
                CenterIV = [1, 0];
            end
            Reward4 = (obj.CenterIV(1)-CenterIV(1));

            reward = Reward1 + Reward2 + Reward3 + Reward4;
            reward = tanh(reward);

            %% Update properties
            obj.ArcSize = length(Archive);
            obj.Zmin    = min(Archive.objs, [], 1);
            obj.Zmax    = max(Archive.objs, [], 1);
            obj.Cmax    = max(max(Archive.cons, 0), [], 1);
            obj.Cmax    = max(obj.Cmax, 1e-6);

            % Normalization
            SConsN1 = max(Archive.cons, 0) ./ repmat(obj.Cmax, obj.ArcSize, 1);
            SConsN1 = sum(SConsN1, 2);
            SObjsN1 = (Archive.objs-repmat(obj.Zmin, obj.ArcSize, 1)) ./ repmat(max(obj.Zmax-obj.Zmin, eps), obj.ArcSize, 1);
            SObjsN1 = sum(SObjsN1, 2);

            obj.Zmin1 = min(SObjsN1);
            obj.Zmax1 = max(SObjsN1);
            obj.Cmax1 = max(SConsN1);

            % Gap solutinon
            CValue     = min(SConsN1);
            OValue     = min(SObjsN1(SConsN1==CValue));
            obj.Center = [OValue, CValue];
           
           %% Centroid point
           [Sobjs, Scons] = obj.Norm(SObjsN1, SConsN1);
            % I Quadrant: +obj, +con
            IndexI = Sobjs >= 0 & Scons > 0;
            if any(IndexI)
                obj.CenterI = mean([Sobjs(IndexI), Scons(IndexI)],1);
            else
                obj.CenterI = [1, 1];
            end
            obj.NRI = sum(IndexI) / obj.ArcSize;

            % II Quadrant: -obj, +con
            IndexII = Sobjs < 0 & Scons > 0;
            if any(IndexII)
                obj.CenterII = mean([Sobjs(IndexII), Scons(IndexII)],1);
            else
                obj.CenterII = [0,0];
            end
            obj.NRII = sum(IndexII) / obj.ArcSize;

            % III Quadrant: -obj, -con
            % empty set

            % IV Quadrant: +obj, -con
            IndexIV = Sobjs > 0 & Scons <= 0;
            if any(IndexIV)
                obj.CenterIV = min([Sobjs(IndexIV), Scons(IndexIV)],[],1);
            else
                obj.CenterIV = [1, 0];
            end
            obj.NRIV = sum(IndexIV) / obj.ArcSize;

            NumPRel = corrcoef(Sobjs,Scons);
            NumPRel = NumPRel(1,2);
            if isnan(NumPRel)
                NumPRel = 0;
            end

            state = [obj.CenterI, obj.CenterII, obj.CenterIV, obj.NRI, obj.NRII, obj.NRIV, NumPRel];
        end


        function [SObjs, SCons] = Norm(obj, Objs, Cons)
            % obj.Center  = [OValue, CValue];
            SObjs = Objs;
            SCons = Cons;
            %% Normalization
            % I Quadrant
            IIndex = Objs>=obj.Center(1) & Cons>obj.Center(2);
            if sum(IIndex)>1
                SObjs(IIndex) = (Objs(IIndex) - obj.Center(1)) / (max(obj.Zmax1 - obj.Center(1), eps));
                SCons(IIndex) = (Cons(IIndex) - obj.Center(2)) / (max(obj.Cmax1 - obj.Center(2), eps));
            end

            % II Quadrant
            IIIndex = Objs<obj.Center(1) & Cons>obj.Center(2);
            if sum(IIIndex)>1
                SObjs(IIIndex) = (Objs(IIIndex) - obj.Zmin1) / (max(obj.Center(1) - obj.Zmin1, eps)) - 1;
                SCons(IIIndex) = (Cons(IIIndex) - obj.Center(2)) / (max(obj.Cmax1 - obj.Center(2), eps));
            end

            % III Quadrant
            IIIIndex = Objs<=obj.Center(1) & Cons<=obj.Center(2);
            IIIIndex(Objs==obj.Center(1) & Cons==obj.Center(2)) = false;
            if sum(IIIIndex)>1
                SObjs(IIIIndex) = (Objs(IIIIndex) - obj.Zmin1) / (max(obj.Center(1) - obj.Zmin1, eps)) - 1;
                SCons(IIIIndex) = (Cons(IIIIndex) + obj.Cmax1) / (max(obj.Center(2) + obj.Cmax1, eps)) - 1;
            end

            % IV Quadrant
            IVIndex = SObjs<obj.Center(1) & Cons<=obj.Center(2);
            if sum(IVIndex)>1
                SObjs(IVIndex) = (Objs(IVIndex) - obj.Center(1)) / (max(obj.Zmax1 - obj.Center(1), eps));
                SCons(IVIndex) = (Cons(IVIndex) + obj.Cmax1) / (max(obj.Center(2) + obj.Cmax1, eps)) - 1;
            end
            if any(SObjs==inf)
                SObjs(SObjs==inf) = 1;
            elseif any(SObjs==-inf)
                SObjs(SObjs==-inf) = 0;
            end

            % Orignal point
            R = find(Objs==obj.Center(1) & Cons==obj.Center(2));
            SCons(R) = 0;
            SObjs(R) = 0;
        end

        function [Population, PopIndex] = FormPop(obj, Archive)
            % Normalization
            SConsN1 = max(Archive.cons, 0) ./ repmat(obj.Cmax, obj.ArcSize, 1);
            SConsN1 = sum(SConsN1, 2);
            SObjsN1 = (Archive.objs-repmat(obj.Zmin, obj.ArcSize, 1)) ./ repmat(obj.Zmax-obj.Zmin, obj.ArcSize, 1);
            SObjsN1 = sum(SObjsN1, 2);
            [Sobjs, Scons] = obj.Norm(SObjsN1, SConsN1);

            % Reverse solutions in II Quadrant
            IndexII = Sobjs < 0 & Scons > 0;
            Sobjs(IndexII) = -Sobjs(IndexII);
            [FrontNo, MaxF] = NDSort([Sobjs, Scons], obj.ArcSize);
            Del = find(FrontNo==MaxF);
            IndexTotal = 1:1:obj.ArcSize;
            PopIndex = RandDel(FrontNo, IndexTotal, obj.ArcSize/2);
            if isscalar(Del)
                PopIndex(end-numel(Del)+1:end) = Del;
            end

            Population = Archive(PopIndex);
        end
    end
end

function PopIndex = RandDel(FrontNo, Index, N)
    Total  = length(FrontNo);
    Num = Total - N;
    if Total >= 2*Num
        Candiates  = randperm(Total, Num*2);
        Candiates1 = Candiates(1:Num);
        Candiates2 = Candiates(Num+1:end);
        Select = [Candiates1(FrontNo(Candiates1)<=FrontNo(Candiates2)),Candiates2(FrontNo(Candiates2)<FrontNo(Candiates1))];
        PopIndex = Index(Select);
        Index(Candiates) = [];
        PopIndex = [PopIndex;Index];
    else
        Candiates  = randperm(Total, N*2);
        Candiates1 = Candiates(1:N);
        Candiates2 = Candiates(N+1:end);
        Select = [Candiates1(FrontNo(Candiates1)<=FrontNo(Candiates2)),Candiates2(FrontNo(Candiates2)<FrontNo(Candiates1))];
        PopIndex = Index(Select);
    end
end