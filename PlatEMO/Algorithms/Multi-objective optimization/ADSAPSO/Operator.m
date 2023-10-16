function Offspring = Operator(Problem,Arc,k,beta,N_a,N_s)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
 
% This function is written by Jianqing Lin
    
    %% Environmental Selection
    [Arc,FrontNo,CrowdDis] = EnvironmentalSelection(Arc,min(N_a,size(Arc,2)));

    %% Parameter setting     
    PopDec = Arc.decs;
    D      = Problem.D;
    M      = Problem.M;
    NP     = Problem.N;

    %% Solution Sort & Selection
    [~,index_FNCD] = sortrows([FrontNo;CrowdDis]',[1,-2]);
    Candidate_D    = PopDec(index_FNCD(1:k),:);
    
    Index_Well = index_FNCD(1:N_s);
    Index_Poor = index_FNCD(end-N_s+1:end);
    
    %% Statistics Mean Value
    Model_Dif      = mean(PopDec(Index_Well,:))-mean(PopDec(Index_Poor,:)); 
    Model_Dif_sort = sort(abs(Model_Dif),'descend');

    %% Selected beta*D Decision Variables
    Index_dif = find(abs(Model_Dif) >= Model_Dif_sort(ceil(beta*D)));
    
    %% Build RBF Models for Low-dimensional Decision Space
    Decs_Surrogate = PopDec(:,Index_dif);
    Objs_Surrogate = Arc.objs;
    
    for i = 1 : M
        RBF_para{i} = RBFCreate(Decs_Surrogate, Objs_Surrogate(:,i), 'gaussian');
    end
   
    %% Reproduction by RBF-assisted PSO
    Population      = EnvironmentalSelection(Arc,NP);
    PopDec          = Population.decs;
    Pop_Surrogate   = Surrogate_individual(PopDec(:,Index_dif),Population.objs);
    Offspring_d     = Reproduction(Problem,Pop_Surrogate,RBF_para,Index_dif);
    
    %% Replacement 
    Offspring_d              = EnvironmentalSelection(Offspring_d,k);
    Candidate_D(:,Index_dif) = Offspring_d.decs;
    Offspring                = Candidate_D;
end