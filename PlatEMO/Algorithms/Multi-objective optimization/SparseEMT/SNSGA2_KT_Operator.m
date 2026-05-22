function [P,Tasks] = SNSGA2_KT_Operator(Problem,Tasks,evaluations)
%% Opeartor of S-NSGA-II
        
        % Generate populations for Knowledge transfer
        P            = [Tasks(1).Pop,Tasks(2).Pop,Tasks(3).Pop];
        [P,FrontNo,CrowdDis] = NSGA2_EnvironmentalSelection(P,Problem.N);
        
        % Knowledge transfer
        maximum = min(Problem.maxFE,Problem.FE + evaluations);
        while Problem.FE < maximum
            ParentPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
            Parents   = P(ParentPool);
            ParentDec1 = [];
            ParentDec2 = [];
            ParentDec3 = [];
            ParentDec4 = [];
            ParentDec5 = [];
            ParentDec6 = [];
            Offspring1 = [];
            Offspring2 = [];
            Offspring3 = [];
            Offspring4 = [];
            Offspring5 = [];
            Offspring6 = [];

            count=1;
            for i=1:2:length(Parents)
                p1=i;
                p2=i+1;
                Dec1 = Parents(p1).dec;
                Dec2 = Parents(p2).dec;
                Dec_high1 = Parents(p1).dec_high;
                Dec_high2 = Parents(p2).dec_high;                

                % Two individuals share the same skill factor
                if Parents(p1).skill_factor==Parents(p2).skill_factor

                    if Parents(p1).skill_factor == 3
                        ParentDec1 = [ParentDec1;Dec1];
                        ParentDec2 = [ParentDec2;Dec2];
                    elseif Parents(p1).skill_factor == 1
                        ParentDec3 = [ParentDec3;Dec1];
                        ParentDec4 = [ParentDec4;Dec2];
                    else
                        ParentDec5 = [ParentDec5;Dec1];
                        ParentDec6 = [ParentDec6;Dec2];
                    end


                else
                    % Two individuals with different skill factors
                    if (Parents(p1).skill_factor == 1 && Parents(p2).skill_factor== 2) || (Parents(p1).skill_factor == 2 && Parents(p2).skill_factor == 1)
                        if rand(1) > 0.5
                            if length(Dec1)>length(Dec2)
                                Dim = Tasks(Parents(p1).skill_factor).Dim;
                                Dec2 = Dec_high2(:,Dim);
                                ParentDec3 = [ParentDec3;Dec1];
                                ParentDec4 = [ParentDec4;Dec2];
                            else
                                Dim = Tasks(Parents(p2).skill_factor).Dim;
                                Dec1 = Dec_high1(:,Dim);
                                ParentDec3 = [ParentDec3;Dec1];
                                ParentDec4 = [ParentDec4;Dec2];
                            end
                        else
                            if length(Dec1)>length(Dec2)
                                Dim = Tasks(Parents(p2).skill_factor).Dim;
                                Dec1 = Dec_high1(:,Dim);
                                ParentDec5 = [ParentDec5;Dec1];
                                ParentDec6 = [ParentDec6;Dec2];
                            else
                                Dim = Tasks(Parents(p1).skill_factor).Dim;
                                Dec2 = Dec_high2(:,Dim);
                                ParentDec5 = [ParentDec5;Dec1];
                                ParentDec6 = [ParentDec6;Dec2];
                            end
                        end
                    end
               


                    if (Parents(p1).skill_factor == 1 && Parents(p2).skill_factor== 3) || (Parents(p1).skill_factor == 3 && Parents(p2).skill_factor == 1)
                        if rand(1) > 0.5
                            if length(Dec1)>length(Dec2)
                                Dim = Tasks(Parents(p1).skill_factor).Dim;
                                Dec2 = Dec_high2(:,Dim);
                                ParentDec1 = [ParentDec1;Dec1];
                                ParentDec2 = [ParentDec2;Dec2];
                            else
                                Dim = Tasks(Parents(p2).skill_factor).Dim;
                                Dec1 = Dec_high1(:,Dim);
                                ParentDec1 = [ParentDec1;Dec1];
                                ParentDec2 = [ParentDec2;Dec2];
                            end
                        else
                            if length(Dec1)>length(Dec2)
                                Dim = Tasks(Parents(p2).skill_factor).Dim;
                                Dec1 = Dec_high1(:,Dim);
                                ParentDec3 = [ParentDec3;Dec1];
                                ParentDec4 = [ParentDec4;Dec2];
                            else
                                Dim = Tasks(Parents(p1).skill_factor).Dim;
                                Dec2 = Dec_high2(:,Dim);
                                ParentDec3 = [ParentDec3;Dec1];
                                ParentDec4 = [ParentDec4;Dec2];
                            end
                        end
                    end




                    if (Parents(p1).skill_factor == 2 && Parents(p2).skill_factor== 3) || (Parents(p1).skill_factor == 3 && Parents(p2).skill_factor == 2)
                        if rand(1) > 0.5
                            if length(Dec1)>length(Dec2)
                                Dim = Tasks(Parents(p1).skill_factor).Dim;
                                Dec2 = Dec_high2(:,Dim);
                                ParentDec1 = [ParentDec1;Dec1];
                                ParentDec2 = [ParentDec2;Dec2];
                            else
                                Dim = Tasks(Parents(p2).skill_factor).Dim;
                                Dec1 = Dec_high1(:,Dim);
                                ParentDec1 = [ParentDec1;Dec1];
                                ParentDec2 = [ParentDec2;Dec2];
                            end
                        else
                            if length(Dec1)>length(Dec2)
                                Dim = Tasks(Parents(p2).skill_factor).Dim;
                                Dec1 = Dec_high1(:,Dim);
                                ParentDec5 = [ParentDec5;Dec1];
                                ParentDec6 = [ParentDec6;Dec2];
                            else
                                Dim = Tasks(Parents(p1).skill_factor).Dim;
                                Dec2 = Dec_high2(:,Dim);
                                ParentDec5 = [ParentDec5;Dec1];
                                ParentDec6 = [ParentDec6;Dec2];
                            end
                        end                        
                    end


                end


                if ~isempty(ParentDec1)   
                    OffDec1       = SNSGA2_SubOperator(Problem,[ParentDec1;ParentDec2],{1,20,1,20,1,20,@spm,@ssbx});
                    OffMask1      = ones(size(OffDec1,1),size(OffDec1,2));
                    OffDec_high1  = Dec_high1;
                    OffDec_high4  = Dec_high2;
                    OffDec_high1(Tasks(3).Dim)  = OffDec1(1,:);
                    OffDec_high4(Tasks(3).Dim)  = OffDec1(2,:);
                    OffMask_high1 = ones(size(OffDec_high1,1),size(OffDec_high1,2));
                    OffMask_high4 = ones(size(OffDec_high4,1),size(OffDec_high4,2));
                    Offspring1    = SOLUTION_SparseEMT(OffDec1(1,:),OffMask1(1,:),OffDec_high1,OffMask_high1,Tasks,3,Problem);
                    Offspring4    = SOLUTION_SparseEMT(OffDec1(2,:),OffMask1(2,:),OffDec_high4,OffMask_high4,Tasks,3,Problem);
                    count = count + 2;
                end

                if ~isempty(ParentDec3)
                    N = size(ParentDec3,1);
                    OffDec2       = SparseEMT_SubOperator(Problem,[ParentDec3;ParentDec4],ones(2*size(ParentDec3,1),size(ParentDec3,2)),Tasks(1).form,Tasks(1).Dim);
                    OffMask2      = ones(size(OffDec2,1),size(OffDec2,2));
                    if N  > 1
                        OffDec_high2  = Dec_high1;
                        OffDec_high5  = Dec_high2;
                        OffDec_high2(Tasks(1).Dim)  = OffDec2(1,:);
                        OffDec_high5(Tasks(1).Dim)  = OffDec2(2,:);
                        OffMask_high2 = ones(size(OffDec_high2,1),size(OffDec_high2,2));
                        OffMask_high5 = ones(size(OffDec_high5,1),size(OffDec_high5,2));
                        Offspring2    = SOLUTION_SparseEMT(OffDec2(1,:),OffMask2(1,:),OffDec_high2,OffMask_high2,Tasks,1,Problem);
                        Offspring5    = SOLUTION_SparseEMT(OffDec2(2,:),OffMask2(2,:),OffDec_high5,OffMask_high5,Tasks,1,Problem);
                        count = count + 2;
                    else                               
                        if rand < 0.5
                            OffDec_high2  = Dec_high1;
                            OffDec_high2(Tasks(1).Dim)  = OffDec2;
                        else
                            OffDec_high2  = Dec_high2;
                            OffDec_high2(Tasks(1).Dim)  = OffDec2;
                        end
                        OffMask_high2 = ones(size(OffDec_high2,1),size(OffDec_high2,2));
                        Offspring2    = SOLUTION_SparseEMT(OffDec2,OffMask2,OffDec_high2,OffMask_high2,Tasks,1,Problem);
                        count = count + 1;
                    end
                end

                if ~isempty(ParentDec5)
                    N = size(ParentDec5,1);
                    OffDec3       = SparseEMT_SubOperator(Problem,[ParentDec5;ParentDec6],ones(2*size(ParentDec5,1),size(ParentDec5,2)),Tasks(2).form,Tasks(2).Dim);
                    OffMask3      = ones(size(OffDec3,1),size(OffDec3,2));
                    if N > 1
                        OffDec_high3  = Dec_high1;
                        OffDec_high6  = Dec_high2;
                        OffDec_high3(Tasks(1).Dim)  = OffDec3(1,:);
                        OffDec_high6(Tasks(1).Dim)  = OffDec3(2,:);
                        OffMask_high3 = ones(size(OffDec_high3,1),size(OffDec_high3,2));
                        OffMask_high6 = ones(size(OffDec_high6,1),size(OffDec_high6,2));
                        Offspring3    = SOLUTION_SparseEMT(OffDec3(1,:),OffMask3(1,:),OffDec_high3,OffMask_high3,Tasks,2,Problem);
                        Offspring6    = SOLUTION_SparseEMT(OffDec3(2,:),OffMask3(2,:),OffDec_high6,OffMask_high6,Tasks,2,Problem);
                        count = count + 2;
                    else
                        if rand < 0.5
                            OffDec_high3  = Dec_high1;
                            OffDec_high3(Tasks(2).Dim)  = OffDec3;
                        else
                            OffDec_high3  = Dec_high2;
                            OffDec_high3(Tasks(2).Dim)  = OffDec3;
                        end
                        OffMask_high3 = ones(size(OffDec_high3,1),size(OffDec_high3,2));
                        Offspring3    = SOLUTION_SparseEMT(OffDec3,OffMask3,OffDec_high3,OffMask_high3,Tasks,2,Problem);
                        count = count + 1;
                    end
                end


                [P,FrontNo,CrowdDis] = NSGA2_EnvironmentalSelection([P,Offspring1,Offspring2,Offspring3,Offspring4,Offspring5,Offspring6],Problem.N);
            end

        end

