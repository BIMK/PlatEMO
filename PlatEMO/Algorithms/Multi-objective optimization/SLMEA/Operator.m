function [OffDec,OffMask,long,numGroup] = Operator(Problem,ParentDec,ParentMask,Fitness,Mix,P,T,Zero,One,Upper,Lower,dt,numGroup,useGPU,GlobalGen)
% The operator of SLMEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Initialization
    [N,D]        = size(ParentDec);
    ClOffMask    = [];
    OrOffMask    = [];
    Parent1Mask  = ParentMask(1:N/2,:);
    Parent2Mask  = ParentMask(N/2+1:end,:);
    Parent11Mask = [];
    Parent12Mask = [];
    Parent21Mask = [];
    Parent22Mask = [];
    
    ClOffDec    = [];
    OrOffDec    = [];
    Parent1Dec  = ParentDec(1:N/2,:);
    Parent2Dec  = ParentDec(N/2+1:end,:);
    Parent11Dec = [];
    Parent12Dec = [];
    Parent21Dec = [];
    Parent22Dec = [];    
    
    long    = 0;
    Fitness = Fitness(Mix);
    
   %% Group
    if GlobalGen >= T
        numGroup = min(200,min(length(Mix),ceil(numGroup*dt)));
        if ~isempty(Mix)
           if ~useGPU
               [index,MAX] = CpuGroup(numGroup,Fitness,length(Mix));
           else
               [index,MAX] = GpuGroup(numGroup,Fitness,length(Mix));
           end
        else
            index = [];
            MAX   = 0;
        end
       if ~useGPU
          Index = zeros(1,length(Lower)); 
       else
          Index = gpuArray.zeros(1,length(Lower)); 
       end
       if ~isempty(Zero)&&~isempty(One)     % Grouping all 1 & all 0
           Index(Zero) = MAX+1;
           Index(One)  = MAX+2;
           MAX = MAX+2;
       elseif isempty(Zero)&&~isempty(One)  % Only Grouping all 1
           Index(One) = MAX+1;
           MAX = MAX+1;
       elseif isempty(One)&&~isempty(Zero)  % Only Grouping all 0
           Index(Zero) = MAX+1;
           MAX = MAX+1;
       end
       Index(Mix) = index;
    end
    %% Group and encoder
    location = rand(1,N/2)<P;
    if ~isempty(find(location>0, 1))&&GlobalGen>=T
       [Parent11Mask]    = Encoder(Parent1Mask(location,:),Index,MAX,useGPU);
       [Parent12Mask]    = Encoder(Parent2Mask(location,:),Index,MAX,useGPU);
       Parent21Mask      = Parent1Mask(~location,:);
       Parent22Mask      = Parent2Mask(~location,:);
       Parent11EncodeDec = (Parent1Dec(location,:)-repmat(Lower,length(find(location==1)),1))./repmat((Upper-Lower),length(find(location==1)),1);
       Parent12EncodeDec = (Parent2Dec(location,:)-repmat(Lower,length(find(location==1)),1))./repmat((Upper-Lower),length(find(location==1)),1);
       Parent11Dec       = EncoderDec(Parent11EncodeDec,Index,MAX);
       Parent12Dec       = EncoderDec(Parent12EncodeDec,Index,MAX);
       Parent21Dec       = Parent1Dec(~location,:);
       Parent22Dec       = Parent2Dec(~location,:);
       if useGPU == 0
           lower = zeros(1,MAX);
           upper = ones(1,MAX);           
       else
           lower = gpuArray.zeros(1,MAX);
           upper = gpuArray.ones(1,MAX);           
       end
    else
        Parent21Mask = Parent1Mask;
        Parent22Mask = Parent2Mask;
        
        Parent21Dec = Parent1Dec;
        Parent22Dec = Parent2Dec;                
    end
    %% Crossover and mutation for mask
    if ~isempty(find(location>0, 1))&&GlobalGen>=T
        ClOffMask = SLMEA_GAhalf([Parent11Mask;Parent12Mask],lower,upper,'binary',useGPU);
        ClOffMask = Decode(ClOffMask,Index);
        long      = length(ClOffMask(:,1));
    end
    OrOffMask = SLMEA_GAhalf([Parent21Mask;Parent22Mask],Lower,Upper,'binary',useGPU);
    OffMask   = [ClOffMask;OrOffMask];
    if any(Problem.encoding~=4)
        if ~isempty(find(location>0, 1))&&GlobalGen>=T 
            ClOffDec = SLMEA_GAhalf([Parent11Dec;Parent12Dec],lower,upper,'real',useGPU);
            ClOffDec = DecodeDec(ClOffDec,Index);
            ClOffDec = ClOffDec.*repmat((Upper-Lower),length(find(location==1)),1)+repmat(Lower,length(find(location==1)),1);
        end
        OrOffDec = SLMEA_GAhalf([Parent21Dec;Parent22Dec],Lower,Upper,'real',useGPU);
        OffDec   = [ClOffDec;OrOffDec];
        OffDec   = min(max(OffDec,repmat(Lower,N/2,1)),repmat(Upper,N/2,1));
        OffDec(:,Problem.encoding==4) = 1;
    else
        if ~useGPU
        	OffDec = ones(N/2,D); 
        else
        	OffDec = gpuArray.ones(N/2,D); 
        end
    end
end

function NewMask= Encoder(Mask,Index,numGroup,useGPU)
    C  = Index;
    Uc = linspace(1,numGroup,numGroup);
    N  = size(Mask,1);
    CC = arrayfun(@(i)find(C==i),Uc,'UniformOutput',false);
    if ~useGPU
        A = cellfun(@(i)mean(Mask(:,i),2)>rand(N,1),CC,'UniformOutput',false);
    else
        A = cellfun(@(i)mean(Mask(:,i),2)>gpuArray.rand(N,1),CC,'UniformOutput',false);
    end
    NewMask = cell2Mat(A);
end

function OMask = Decode(NewMask,Index)
	OMask = NewMask(:,Index);
end
function NewDec = EncoderDec(Dec,Index,numGroup)
    C      = Index;
    Uc     = linspace(1,numGroup,numGroup);
    CC     = arrayfun(@(i)find(C==i),Uc,'UniformOutput',false);
    A      = cellfun(@(i)mean(Dec(:,i),2),CC,'UniformOutput',false);
    NewDec = cell2Mat(A);
end

function ODec = DecodeDec(NewDec,Index)
	ODec = NewDec(:,Index);
end

function m = cell2Mat(c)
    rows = size(c,1);
    cols = size(c,2);   
    if rows < cols
        m = cell(rows,1);
        for n = 1 : rows
            m{n} = cat(2,c{n,:});
        end
        m = cat(1,m{:});
    else
        m = cell(1, cols);
        for n = 1 : cols
            m{n} = cat(1,c{:,n});
        end    
        m = cat(2,m{:});
    end
end