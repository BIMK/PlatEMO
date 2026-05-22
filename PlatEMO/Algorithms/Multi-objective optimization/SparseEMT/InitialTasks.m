function [Fitness,FitnessObj,Tasks,Dec,Mask] = InitialTasks(Problem,no_of_tasks,solver)
%% Update multiple tasks

     % Calculate the fitness of each decision variable
     [Fitness,FitnessObj,TDec,TMask,Dim] = FitnessInit(Problem,solver);

     % Classify decision variables
     REAL          = any(Problem.encoding==1);
     idx_union     = any(Dim,1);
     idx_intersect = all(Dim,1);
     
     
     % Generate initial population
     if solver ~= 2
         Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
         Dec(:,Problem.encoding==4) = 1;
         Mask = false(Problem.N,Problem.D);
         for i = 1 : Problem.N
             Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
         end
         if REAL && Problem.D >= 2000
             AllSample   = randperm(length(TMask));
             FinalSample = AllSample(1:Problem.D);
             TDec        = TDec(FinalSample,:);
             TMask       = TMask(FinalSample,:);
         end
         Dec = [Dec;TDec];
         Mask = [Mask;TMask];
     else
         sampling_method = {@SNSGA2_vssps, 0.75, 1};
         sampler = sampling_method{1};
         lowerBound = sampling_method{2};
         upperBound = sampling_method{3};
         Dec =  sampler(Problem, lowerBound, upperBound);
         Mask = ones(size(Dec,1),size(Dec,2));
     end

     % Initialize or update tasks
     for i = 1:no_of_tasks
         index = [];
         switch i
         case 1
             index = find(idx_union==1);
             index = index';
             Tasks(i).form = 3;
         case 2
             index = find(idx_intersect==1);
             index = index';
             Tasks(i).form = 2;
         case no_of_tasks
             index = 1:1:Problem.D;
             index = index';
             Tasks(i).form = 1;
         end
         Tasks(i).skill_factor = i;
         Tasks(i).Dim          = index;
         if REAL
             Tasks(i).upper        = Problem.upper(index);
             Tasks(i).lower        = Problem.lower(index);
         else
             Tasks(i).upper        = ones(1,length(index));
             Tasks(i).lower        = zeros(1,length(index));
         end
         Tasks(i).Pop          = [];
     end
end


function [Fitness,FitnessObj,TDec,TMask,Dim] = FitnessInit(Problem,solver)
%% Calculate the fitness of each decision variable and generate initial population

      % Calculate the fitness of each decision variable
      TDec    = [];
      TMask   = [];
      REAL    = any(Problem.encoding==1);
      if REAL       
          FitnessObj    = zeros(5,Problem.D);
          Dim           = zeros(5,Problem.D);
          Mask_Temp     = zeros(1,Problem.D);
          Interval = (Problem.upper-Problem.lower)./5;
          for i = 1 : 5  
              if solver~= 3
                  Dec  = unifrnd(repmat(Problem.lower + Interval*(i-1),Problem.D+1,1),repmat(Problem.lower + Interval*(i),Problem.D+1,1));
                  Mask = [Mask_Temp;eye(Problem.D)];
                  Population   = Problem.Evaluation(Dec.*Mask);
                  TDec         = [TDec;Dec];
                  TMask        = [TMask;Mask];
                  FrontNo      = NDSort([Population.objs,Population.cons],inf);
                  FitnessObj (i,:) = FitnessObj(i,:) + FrontNo(2:end);
                  Dim(i,:) = FrontNo(1)>=FrontNo(2:end);
              else
                  for j = 1 : 2
                      Dec  = unifrnd(repmat(Problem.lower + Interval*(i-1),Problem.D+1,1),repmat(Problem.lower + Interval*(i),Problem.D+1,1));
                      Mask = [Mask_Temp;eye(Problem.D)];
                      Population   = Problem.Evaluation(Dec.*Mask);
                      TDec         = [TDec;Dec];
                      TMask        = [TMask;Mask];
                      FrontNo      = NDSort([Population.objs,Population.cons],inf);
                      FitnessObj (i,:) = FitnessObj(i,:) + FrontNo(2:end);
                      Dim(i,:) = FrontNo(1)>=FrontNo(2:end);
                  end
              end
          end
          Fitness = sum(FitnessObj);  
      else
          FitnessObj = zeros(1,Problem.D);
          Dim        = zeros(1,Problem.D);
          Mask_Temp  = zeros(1,Problem.D);
          for i = 1 
              Dec          = ones(Problem.D+1,Problem.D);
              Mask         = [Mask_Temp;eye(Problem.D)];
              Population   = Problem.Evaluation(Dec.*Mask);
              TDec         = [TDec;Dec];
              TMask        = [TMask;Mask];
              FrontNo      = NDSort([Population.objs,Population.cons],inf);
              FitnessObj (i,:) = FitnessObj(i,:) + FrontNo(2:end); 
              Dim(i,:) = FrontNo(1)>=FrontNo(2:end);
          end
          Fitness = FitnessObj;
      end
end