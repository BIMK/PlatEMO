function Tasks = DividePop(Problem,Tasks,RefPop,solver,Parameter)
%% Opeartor of SparseEMT    
        REAL        = any(Problem.encoding==1);
        N_auxiliary = floor(Problem.N/4);
        switch solver
            case 1
                RefPop1 = MSKEA_EnvironmentalSelection(RefPop,N_auxiliary);
            case 2
                RefPop1 = NSGA2_EnvironmentalSelection(RefPop,N_auxiliary);
            case 3
                RefPop1 = MGCEA_EnvironmentalSelection(RefPop,N_auxiliary);
        end  
        % Generate population for the auxiliary Task1
        if REAL
            Mask_high     = RefPop1.masks_high;
            Mask          = Mask_high(:,Tasks(1).Dim);
            Dec_high      = RefPop1.decss_high;
            Dec           = Dec_high(:,Tasks(1).Dim);
        else
            Mask_high     = RefPop1.masks_high;
            Mask          = Mask_high(:,Tasks(1).Dim);
            Dec_high      = ones(size(Mask,1),size(Mask_high,2));
            Dec           = ones(size(Mask,1),length(Tasks(1).Dim));
        end
        Tasks(1).Pop  = SOLUTION_SparseEMT(Dec,Mask,Dec_high,Mask_high,Tasks,1,Problem);
       

        % Generate population for the auxiliary Task2
        if REAL
            Mask_high     = RefPop1.masks_high; 
            Mask_high(:,Tasks(2).Dim) = 1;
            Mask          = ones(size(Mask_high,1),length(Tasks(2).Dim));
            Dec_high      = RefPop1.decss_high;
            Dec           = Dec_high(:,Tasks(2).Dim);
            Tasks(2).Pop  = SOLUTION_SparseEMT(Dec,Mask,Dec_high,Mask_high,Tasks,2,Problem);
        else
            Tasks(2).Pop = [];
        end
       
        
        % Generate population for the original Task
        N_original  = Problem.N;
        Tasks(3).Pop  = RefPop(RefPop.skill_factors == 3);
        if length(Tasks(3).Pop) < N_original
            inx = find(RefPop.skill_factors~=3);
            k = min(N_original-length(Tasks(3).Pop),length(inx));
            inx=inx(randperm(length(inx)));
            Dec           = RefPop(inx(1:k)).decss_high;
            Mask          = RefPop(inx(1:k)).masks_high;
            P = [Tasks(3).Pop,SOLUTION_SparseEMT(Dec,Mask,Dec,Mask,Tasks,3,Problem)];
            maximum = Problem.FE + 5e-2*Problem.maxFE;
            switch solver
                case 1
                    [sv,pv,Last_temp_num] = deal(Parameter{:});
                    [P,FrontNo,CrowdDis] = MSKEA_EnvironmentalSelection(P,N_original);
                    while Problem.FE < maximum
                        ParentPool  = TournamentSelection(2,N_original,FrontNo,-CrowdDis);
                        Parents     = P(ParentPool);
                        %-------------update fv-----------%
                        delta = Problem.FE/Problem.maxFE;
                        if delta<0.618
                            fv = std(P(FrontNo==1).decss_high.*P(FrontNo==1).masks_high,0,1);
                            temp = P.masks_high;
                            fv(:,Problem.encoding==4) = sum(temp(FrontNo==1,Problem.encoding==4),1);
                        end
                        %-------------update sv-----------%
                        First_Mask=P(FrontNo==1).masks_high;
                        [temp_num,~]=size(First_Mask);
                        temp_vote=sum(First_Mask,1);
                        sv(1,:)=(Last_temp_num/(Last_temp_num+temp_num))*sv(1,:)+(temp_num/(Last_temp_num+temp_num))*(temp_vote/temp_num);
                        Last_temp_num=temp_num;
                        %-------------update pv by sv-----------%
                        if delta<0.618
                            pv=pv.*(1-sv)*sqrt(delta)+pv;
                        end
                        delta = Problem.FE/Problem.maxFE;
                        if  (delta/0.618) < 0.618
                            [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,Parents.decss,Parents.masks,pv,fv,delta);
                        elseif (delta/0.618)>=0.618 && delta< 0.618
                            if rand < 0.5
                                [OffDec,OffMask] = MSKEA_Operator_sv(Problem,Parents.decss,Parents.masks,sv);
                            else
                                [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,Parents.decss,Parents.masks,pv,fv,delta);
                            end
                        else
                            [OffDec,OffMask] = MSKEA_Operator_sv(Problem,Parents.decss,Parents.masks,sv);
                        end
                        Offspring   = SOLUTION_SparseEMT(OffDec,OffMask,OffDec,OffMask,Tasks,3,Problem);
                        [P,FrontNo,CrowdDis] = MSKEA_EnvironmentalSelection([P,Offspring],N_original);
                    end
                case 2
                    [P,FrontNo,CrowdDis] = NSGA2_EnvironmentalSelection(P,N_original);
                    while Problem.FE < maximum
                        ParentPool = TournamentSelection(2,N_original,FrontNo,-CrowdDis);
                        Parents    = P(ParentPool);
                        OffDec     = SNSGA2_SubOperator(Problem,Parents.decs,{1,20,1,20,1,20,@spm,@ssbx});
                        OffMask    = ones(size(OffDec,1),size(OffDec,2));
                        Offspring  = SOLUTION_SparseEMT(OffDec,OffMask,OffDec,OffMask,Tasks,3,Problem);
                        [P,FrontNo,CrowdDis] = NSGA2_EnvironmentalSelection([P,Offspring],N_original);
                    end
                case 3
                    [Fitness,SparseRate,NearStage,FitnessLayer,LayerMax] = deal(Parameter{:});
                    [P,Spea2Fit] = MGCEA_EnvironmentalSelection(P,N_original);
                    while Problem.FE < maximum
                        ParentPool = TournamentSelection(2,N_original,Spea2Fit);
                        Parents    = P(ParentPool);
                        [NearStage,Fitness,FitnessLayer,LayerMax] = MGCEA_ControlStage(SparseRate,NearStage,P.masks_high,Fitness,FitnessLayer,LayerMax,Problem);
                        [OffDec,OffMask] = MGCEA_SubOperator(Problem,Parents.decss,Parents.masks,FitnessLayer,LayerMax);
                        Offspring  = SOLUTION_SparseEMT(OffDec,OffMask,OffDec,OffMask,Tasks,3,Problem);
                        [P,Spea2Fit] = MGCEA_EnvironmentalSelection([P,Offspring],N_original);
                    end
            end
            Tasks(3).Pop = P;
        end
end

