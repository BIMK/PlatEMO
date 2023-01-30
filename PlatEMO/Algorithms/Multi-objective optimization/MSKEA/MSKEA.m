classdef MSKEA < ALGORITHM
 % <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
 % Multi-stage knowledge-guided evolutionary algorithm
 
%------------------------------- Reference --------------------------------
% Z. Ding, L. Chen, D. Sun, and X. Zhang, A multi-stage knowledge-guided
% evolutionary algorithm for sparse multi-objective optimization problems,
% Swarm and Evolutionary Computation, 2022, 73: 101119.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lei Chen

    methods
         function main(Algorithm,Problem)
            %% pv and population initialization  
            TDec    = [];
            TMask   = [];
            TempPop = [];
            pv      = zeros(1,Problem.D);
            for i = 1 : 1+4*any(Problem.encoding~=4)
                Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                Dec(:,Problem.encoding==4) = 1;
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                pv         = pv + NDSort([Population.objs,Population.cons],inf);
            end
            Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Dec(:,Problem.encoding==4) = 1;
            Mask = false(Problem.N,Problem.D);
            
            %% Mask initialization by pv
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),pv)) = 1;
            end
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);
            sv=zeros(1,Problem.D);
            Last_temp_num=0;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                %-------------update fv-----------%
                delta= Problem.FE/Problem.maxFE;
                if delta<0.618
                    fv = std(Population(FrontNo==1).decs,0,1);
                    fv(:,Problem.encoding==4) = sum(Mask(FrontNo==1,Problem.encoding==4),1);
                end
                %-------------update sv-----------%
                First_Mask=Mask(FrontNo==1,:);
                [temp_num,~]=size(First_Mask);
                temp_vote=sum(First_Mask,1);
                sv(1,:)=(Last_temp_num/(Last_temp_num+temp_num))*sv(1,:)+(temp_num/(Last_temp_num+temp_num))*(temp_vote/temp_num);
                Last_temp_num=temp_num;
                %-------------update pv by sv-----------%
                if delta<0.618
                    pv=pv.*(1-sv)*sqrt(delta)+pv;
                end
                %--------------------------------------%
                if  (delta/0.618) < 0.618
                    [OffDec,OffMask] = Operator_pvfv(Problem,Dec(MatingPool,:),Mask(MatingPool,:),pv,fv,delta);
                elseif (delta/0.618)>=0.618 && delta< 0.618
                    if rand < 0.5
                        [OffDec,OffMask] = Operator_sv(Problem,Dec(MatingPool,:),Mask(MatingPool,:),sv);
                    else
                        [OffDec,OffMask] = Operator_pvfv(Problem,Dec(MatingPool,:),Mask(MatingPool,:),pv,fv,delta);
                    end
                else
                    [OffDec,OffMask] = Operator_sv(Problem,Dec(MatingPool,:),Mask(MatingPool,:),sv);
                end
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
            end
        end
    end
end