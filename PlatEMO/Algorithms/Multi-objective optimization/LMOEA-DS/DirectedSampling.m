function GuidingSolution = DirectedSampling(Problem,Population,Ns,Nw,RefV)
% Acquiring Guiding Solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Shufen Qin
% E-mail: shufen.qin@stu.tyust.edu.cn

    %% Classter the reference vectors
    BoundRefV               = eye(Problem.M);
    BoundRefV(BoundRefV==0) = 10e-7;
    [~,CenterRefV,~,~]      = kmeans(RefV,Nw);
    DirectRefV              = [BoundRefV;CenterRefV];
    Nw                      = size(DirectRefV,1);
    
    %% Identify guiding directions
    Best       = GenerateRepresetativeSolution(Population.objs,DirectRefV);
    PopDec     = Population.decs;
    BestX      = PopDec(Best,:);
    Upper      = Problem.upper; 
    Lower      = Problem.lower;
    Directnorm = [sqrt(sum((BestX - repmat(Lower,Nw,1)).^2,2));sqrt(sum((BestX - repmat(Upper,Nw,1)).^2,2))];
    Direction  = [BestX - repmat(Lower,Nw,1);BestX - repmat(Upper,Nw,1)]./repmat(Directnorm,1,Problem.D);
    
    %% Generate guiding solutions
    Intervalmax     = sqrt(sum((Upper-Lower).^2,2));
    Intervalmin       = 0;
    Nw              = 2*Nw;
    RandSample      = Intervalmin + rand(Ns,Nw)*(Intervalmax-Intervalmin);
    SampleSolution  = GenerateSampleSolution(Problem,RandSample,Direction);
    GuidingSolution = SampleSolution((NDSort(SampleSolution.objs,1)==1));
    
end

function Best = GenerateRepresetativeSolution(Obj,RefV)
% Find out respective solutions

    %% Normalization
    np = size(Obj,1);
    Obj = (Obj-repmat(min(Obj),np,1))./(repmat(max(Obj),np,1)-repmat(min(Obj),np,1));
    Nr   = size(RefV,1);
    Best = zeros(Nr,1);

    %% Assign individuals for each reference vector
    Cosine        = 1-pdist2(Obj,RefV,'cosine');
    [~,associate] = max(Cosine,[],2);

    Indflag = zeros(np,1);
    current = cell(Nr,1);
    for i = 1 : Nr
        current{i,1} = find(associate == i);
    end  
    for i = 1 : Nr
        if length(current{i,1})>1
            normf = sqrt(sum(Obj(current{i,1},:).^2,2));
            normRefV = sqrt(sum(repmat(RefV(i,:),length(current{i,1}),1).^2,2));
            CosinefRefV = sum(Obj(current{i,1},:).*repmat(RefV(i,:),length(current{i,1}),1),2)./normRefV./normf;
            d1 = normf .* CosinefRefV;
            [~,ind] = sort(d1,'ascend');
            Best(i,1) = current{i,1}(ind(1));
            Indflag(current{i,1}(ind(1)),1) = 1;
        elseif length(current{i,1})==1
            Best(i,1) = current{i,1}(1);
            Indflag(current{i,1}(1),1) = 1;
        end
    end
    for i = 1:Nr
        if isempty(current{i,1})
            [~,indCon] = sort(Cosine(:,i),'descend');
            k = 1;
            if length(indCon) > Nr
                while Indflag(indCon(k),1) == 1
                    k=k+1;
                end
                Best(i,1) = indCon(k);
                Indflag(indCon(k),1) = 1;
            else
                Best(i,1) = indCon(1);
            end
        end
    end
end

function SampleSolution = GenerateSampleSolution(Problem,RandSample,Direct)
% Generate some sample solutions along with the guiding directions

   [Ns,Nw] = size(RandSample);
   Nw = Nw/2;
   SampleSolution = [];
   for i = 1:Ns
       PopX = [repmat(Problem.lower,Nw,1) + repmat(RandSample(i,1:Nw)',1,Problem.D).* Direct(1:Nw,:);...
           repmat(Problem.upper,Nw,1) + repmat(RandSample(i,Nw+1:end)',1,Problem.D).* Direct(Nw+1:end,:)];
       PopX = max(min(repmat(Problem.upper,size(PopX,1),1),PopX),repmat(Problem.lower,size(PopX,1),1));
       SampleSolutiontemp = Problem.Evaluation(PopX);
       SampleSolution = [SampleSolution,SampleSolutiontemp];
   end
end