function [P,Tasks,sv,pv,Last_temp_num] = MSKEA_KT_Operator(Problem,Tasks,evaluations,rmp,sv,pv,Last_temp_num)
%% Opeartor of MSKEA
        
        % Generate populations for Knowledge transfer
        P  = [Tasks(1).Pop,Tasks(2).Pop,Tasks(3).Pop];
        [P,FrontNo,CrowdDis] = MSKEA_EnvironmentalSelection(P,Problem.N);
        
        % Knowledge transfer
        maximum = min(Problem.maxFE,Problem.FE + evaluations);
        while Problem.FE < maximum
            ParentPool  = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
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

            count=1;
            for i=1:2:length(Parents)
                p1=i;
                p2=i+1;
                Dec1 = Parents(p1).dec;
                Dec2 = Parents(p2).dec;
                Dec_high1 = Parents(p1).dec_high;
                Dec_high2 = Parents(p2).dec_high;
                Mask1 = Parents(p1).mask;
                Mask2 = Parents(p2).mask;
                Mask_high1 = Parents(p1).mask_high;
                Mask_high2 = Parents(p2).mask_high;
                % Two individuals share the same skill factor
                if Parents(p1).skill_factor==Parents(p2).skill_factor
                    if Parents(p1).skill_factor == 3
                        if  (delta/0.618) < 0.618
                            [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,[Dec1;Dec2],[Mask1;Mask2],pv,fv,delta);
                        elseif (delta/0.618)>=0.618 && delta< 0.618
                            if rand < 0.5
                                [OffDec,OffMask] = MSKEA_Operator_sv(Problem,[Dec1;Dec2],[Mask1;Mask2],sv);
                            else
                                [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,[Dec1;Dec2],[Mask1;Mask2],pv,fv,delta);
                            end
                        else
                            [OffDec,OffMask] = MSKEA_Operator_sv(Problem,[Dec1;Dec2],[Mask1;Mask2],sv);
                        end
                    else
                        [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p1).skill_factor).form,Tasks(Parents(p1).skill_factor).Dim);
                    end
                    skill_factor= Parents(p1).skill_factor;
                else
                    % Two individuals with different skill factors
                    if rand(1)<rmp
                        if (Parents(p1).skill_factor == 1 && Parents(p2).skill_factor== 2) || (Parents(p1).skill_factor == 2 && Parents(p2).skill_factor == 1)
                            if rand(1) > 0.5 
                                if length(Dec1) > length(Dec2)
                                    Dim = Tasks(Parents(p1).skill_factor).Dim;
                                    Dec2 = Dec_high2(:,Dim);
                                    Mask2 = Mask_high2(:,Dim);
                                    [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p1).skill_factor).form,Tasks(Parents(p1).skill_factor).Dim);
                                    skill_factor = Parents(p1).skill_factor;
                                else
                                    Dim = Tasks(Parents(p2).skill_factor).Dim;
                                    Dec1 = Dec_high1(:,Dim);
                                    Mask1 = Mask_high1(:,Dim);
                                    [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p2).skill_factor).form,Tasks(Parents(p2).skill_factor).Dim);
                                    skill_factor = Parents(p2).skill_factor;
                                end
                            else
                                if length(Dec1) > length(Dec2)
                                    Dim = Tasks(Parents(p2).skill_factor).Dim;
                                    Dec1 = Dec_high1(:,Dim);
                                    Mask1 = Mask_high1(:,Dim);
                                    [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p2).skill_factor).form,Tasks(Parents(p2).skill_factor).Dim);
                                    skill_factor = Parents(p2).skill_factor;
                                else
                                    Dim = Tasks(Parents(p1).skill_factor).Dim;
                                    Dec2 = Dec_high2(:,Dim);
                                    Mask2 = Mask_high2(:,Dim);
                                    [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p1).skill_factor).form,Tasks(Parents(p1).skill_factor).Dim);
                                    skill_factor = Parents(p1).skill_factor;
                                end                                
                            end
                        end

                        if (Parents(p1).skill_factor == 1 && Parents(p2).skill_factor== 3) || (Parents(p1).skill_factor == 3 && Parents(p2).skill_factor == 1) ...
                            || (Parents(p1).skill_factor == 2 && Parents(p2).skill_factor== 3) || (Parents(p1).skill_factor == 3 && Parents(p2).skill_factor == 2)
                            if rand(1)<0.5
                                if length(Dec1)>length(Dec2)
                                    Dim = Tasks(Parents(p1).skill_factor).Dim;
                                    Dec2 = Dec_high2(:,Dim);
                                    Mask2 = Mask_high2(:,Dim);
                                    if  (delta/0.618) < 0.618
                                        [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,[Dec1;Dec2],[Mask1;Mask2],pv,fv,delta);
                                    elseif (delta/0.618)>=0.618 && delta< 0.618
                                        if rand < 0.5
                                            [OffDec,OffMask] = MSKEA_Operator_sv(Problem,[Dec1;Dec2],[Mask1;Mask2],sv);
                                        else
                                            [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,[Dec1;Dec2],[Mask1;Mask2],pv,fv,delta);
                                        end
                                    else
                                        [OffDec,OffMask] = MSKEA_Operator_sv(Problem,[Dec1;Dec2],[Mask1;Mask2],sv);
                                    end
                                    skill_factor= Parents(p1).skill_factor;
                                else
                                    Dim = Tasks(Parents(p2).skill_factor).Dim;
                                    Dec1 = Dec_high1(:,Dim);
                                    Mask1 = Mask_high1(:,Dim);
                                    if  (delta/0.618) < 0.618
                                        [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,[Dec1;Dec2],[Mask1;Mask2],pv,fv,delta);
                                    elseif (delta/0.618)>=0.618 && delta< 0.618
                                        if rand < 0.5
                                            [OffDec,OffMask] = MSKEA_Operator_sv(Problem,[Dec1;Dec2],[Mask1;Mask2],sv);
                                        else
                                            [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,[Dec1;Dec2],[Mask1;Mask2],pv,fv,delta);
                                        end
                                    else
                                        [OffDec,OffMask] = MSKEA_Operator_sv(Problem,[Dec1;Dec2],[Mask1;Mask2],sv);
                                    end
                                    skill_factor= Parents(p2).skill_factor;
                                end
                  
                            else
                                if length(Dec1) > length(Dec2)
                                    Dim = Tasks(Parents(p2).skill_factor).Dim;
                                    Dec1 = Dec_high1(:,Dim);
                                    Mask1 = Mask_high1(:,Dim);
                                    [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p2).skill_factor).form,Tasks(Parents(p2).skill_factor).Dim);
                                    skill_factor= Parents(p2).skill_factor;
                                else
                                    Dim = Tasks(Parents(p1).skill_factor).Dim;
                                    Dec2 = Dec_high2(:,Dim);
                                    Mask2 = Mask_high2(:,Dim);
                                    [OffDec,OffMask] = SparseEMT_SubOperator(Problem,[Dec1;Dec2],[Mask1;Mask2],Tasks(Parents(p1).skill_factor).form,Tasks(Parents(p1).skill_factor).Dim);
                                    skill_factor= Parents(p1).skill_factor;
                                end
                            end
                        end

                                

                   
                           
                    else
                        if (Parents(p1).skill_factor == 1 && Parents(p2).skill_factor== 2) || (Parents(p1).skill_factor == 2 && Parents(p2).skill_factor == 1)
                            [OffDec1,OffMask1] = SparseEMT_SubOperator(Problem,Dec1,Mask1,Tasks(Parents(p1).skill_factor).form,Tasks(Parents(p1).skill_factor).Dim);
                            [OffDec2,OffMask2] = SparseEMT_SubOperator(Problem,Dec2,Mask2,Tasks(Parents(p2).skill_factor).form,Tasks(Parents(p2).skill_factor).Dim);
                        end


                        if (Parents(p1).skill_factor == 1 && Parents(p2).skill_factor== 3) || (Parents(p1).skill_factor == 3 && Parents(p2).skill_factor == 1)...
                           || (Parents(p1).skill_factor == 2 && Parents(p2).skill_factor== 3) || (Parents(p1).skill_factor == 3 && Parents(p2).skill_factor == 2)                
                            if length(Dec1)>length(Dec2)
                                if  (delta/0.618) < 0.618
                                    [OffDec1,OffMask1] = MSKEA_Operator_pvfv(Problem,Dec1,Mask1,pv,fv,delta);
                                elseif (delta/0.618)>=0.618 && delta< 0.618
                                    if rand < 0.5
                                        [OffDec1,OffMask1] = MSKEA_Operator_sv(Problem,Dec1,Mask1,sv);
                                    else
                                        [OffDec1,OffMask1] = MSKEA_Operator_pvfv(Problem,Dec1,Mask1,pv,fv,delta);
                                    end
                                else
                                    [OffDec1,OffMask1] = MSKEA_Operator_sv(Problem,Dec1,Mask1,sv);
                                end
                                [OffDec2,OffMask2] = SparseEMT_SubOperator(Problem,Dec2,Mask2,Tasks(Parents(p2).skill_factor).form,Tasks(Parents(p2).skill_factor).Dim);
                            else
                                [OffDec1,OffMask1] = SparseEMT_SubOperator(Problem,Dec1,Mask1,Tasks(Parents(p1).skill_factor).form,Tasks(Parents(p1).skill_factor).Dim);
                                if  (delta/0.618) < 0.618
                                    [OffDec2,OffMask2] = MSKEA_Operator_pvfv(Problem,Dec2,Mask2,pv,fv,delta);
                                elseif (delta/0.618)>=0.618 && delta< 0.618
                                    if rand < 0.5
                                        [OffDec2,OffMask2] = MSKEA_Operator_sv(Problem,Dec2,Mask2,sv);
                                    else
                                        [OffDec2,OffMask2] = MSKEA_Operator_pvfv(Problem,Dec2,Mask2,pv,fv,delta);
                                    end
                                else
                                    [OffDec2,OffMask2] = MSKEA_Operator_sv(Problem,Dec2,Mask2,sv);
                                end
                            end
                        end

                        skill_factor= [Parents(p1).skill_factor,Parents(p2).skill_factor];
                    end
                end
                if length(skill_factor) == 2
                    OffDec_high1  = Dec_high1;
                    OffDec_high2  = Dec_high2;
                    OffMask_high1 = Mask_high1;
                    OffMask_high2 = Mask_high2;
                    OffDec_high1(Tasks(skill_factor(1)).Dim)  = OffDec1;
                    OffDec_high2(Tasks(skill_factor(2)).Dim)  = OffDec2;
                    OffMask_high1(Tasks(skill_factor(1)).Dim)  = OffMask1;
                    OffMask_high2(Tasks(skill_factor(2)).Dim)  = OffMask2;
                    Offspring(count)   = SOLUTION_SparseEMT(OffDec1,OffMask1,OffDec_high1,OffMask_high1,Tasks,skill_factor(1),Problem);
                    Offspring(count+1) = SOLUTION_SparseEMT(OffDec2,OffMask2,OffDec_high2,OffMask_high2,Tasks,skill_factor(2),Problem);
                    count=count+2;
                else
                    if rand < 0.5
                        OffDec_high  = Dec_high1;
                        OffMask_high = Mask_high1;
                        OffDec_high(Tasks(skill_factor).Dim)  = OffDec;
                        OffMask_high(Tasks(skill_factor).Dim)  = OffMask;
                    else
                        OffDec_high  = Dec_high2;
                        OffMask_high = Mask_high2;
                        OffDec_high(Tasks(skill_factor).Dim)  = OffDec;
                        OffMask_high(Tasks(skill_factor).Dim)  = OffMask;
                    end
                    Offspring(count)   = SOLUTION_SparseEMT(OffDec,OffMask,OffDec_high,OffMask_high,Tasks,skill_factor,Problem);
                    count=count+1;
                end
            end
            [P,FrontNo,CrowdDis] = MSKEA_EnvironmentalSelection([P,Offspring],Problem.N);
        end
end

