function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,FitnessLayer,LayerMax)
% The operator of SparseEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [N,D]       = size(ParentDec);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);

    OffMask        = Parent1Mask;
    SparseRate1    = sum(Parent1Mask,2);
    SparseRate2    = sum(Parent2Mask,2);
    Rate           = SparseRate1./(SparseRate1 + SparseRate2);
    index          = rand(N/2,D) > repmat(Rate/2,1,D);
    OffMask(index) = Parent2Mask(index);
    
    for i = 1 : N/2
        PointUp   = 1;
        PointDown = LayerMax;
        for j = 1 : LayerMax      
            TargetUpLayer   = find(FitnessLayer == PointUp);
            TargetUp        = TargetUpLayer(OffMask(i,TargetUpLayer) == 0);
            TargetDownLayer = find(FitnessLayer == PointDown);           
            TargetDown      = TargetDownLayer(OffMask(i,TargetDownLayer) == 1);
            if rand < 0.5
                if ~isempty(TargetUp)
                    if rand < 0.5
                        TargetUp = datasample(TargetUp,ceil(length(TargetUp)/2));
                        OffMask(i,TargetUp) = 1;
                    end
                end
                if rand < 0.5
                    PointUp = PointUp + 1;
                else
                    break;
                end
            else
                if ~isempty(TargetDown)
                    if rand < 0.5
                        TargetDown = datasample(TargetDown,ceil(length(TargetDown)/2));
                        OffMask(i,TargetDown) = 0;
                    end
                end
                if rand < 0.5
                    PointDown = PointDown - 1;
                else
                    break;
                end
            end
            if PointUp >= PointDown
                break;
            end            
        end
    end

    %% Crossover and mutation for dec
    if any(Problem.encoding~=4)
        OffDec = OperatorGAhalf(Problem,ParentDec);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(N/2,D);
    end
end