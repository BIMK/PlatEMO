function MatingPool = MatingSelection(Population,ArcPop,N)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    MatingPool = [];
    if length(ArcPop) < N
        SelectedIndex = TournamentSelection(2,N,-sum(max(0,Population.cons),2));
        MatingPool    = Population(SelectedIndex);
    else
        AllPop = [ Population,ArcPop]; 
        Zmin       = min(AllPop.objs,[],1);
        PopObj = (AllPop.objs-repmat(Zmin,length(AllPop.objs),1))./(repmat(max(AllPop.objs),length(AllPop.objs),1)-repmat(Zmin,length(AllPop.objs),1)+1e-10)+1e-10;
        Cosine   = 1 - pdist2(PopObj,PopObj,'cosine');
        Cosine   = Cosine.*(1-eye(size(PopObj,1)));

        Temp     = sort(-Cosine,2);
        [~,Rank] = sortrows(Temp);

        CV1 = sum(max(0,Population.cons),2);
        CV2 = sum(max(0,ArcPop.cons),2);

        Angle1 = Rank(1:N);
        Angle2 = Rank(N+1:length(AllPop));

        Index1 = randi(N,1,N);
        Index2 = randi(length(ArcPop),1,N);

        i = 0;
        while length(MatingPool)< N  
            if CV1(Index1(i+1))< CV2(Index2(i+1))     
                MatingPool = [MatingPool,Population(Index1(i+1))];
            else
                MatingPool = [MatingPool,ArcPop(Index2(i+1))];
            end
            if Angle1(Index1(i+2))< Angle2(Index2(i+2))
                MatingPool = [MatingPool,Population(Index1(i+2))];
            else
                MatingPool = [MatingPool,ArcPop(Index2(i+2))];
            end    
            i = i + 2 ;
        end
    end
end