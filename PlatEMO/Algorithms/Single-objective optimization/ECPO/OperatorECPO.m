function ECP_NEW = OperatorECPO(Problem,Strategy,ECP,Archive,pop_fac,nECPI)
%OperatorECPO - The operator of Electric Charged Particles Optimization.
%
%   ECP_NEW = OperatorECPO(Problem,Strategy,ECP,Archive,pop_fac,nECPI) 
%
%   Example:
%       Off = OperatorECPO(Problem,Strategy,ECP,Archive,pop_fac,nECPI)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by H. R. E. H. Bouchekara (email: bouchekara.houssem@gmail.com)

    %% Parameter setting
    ECPDec             = ECP.decs;
    ArchiveDEC         = Archive.decs;
    [nECP,PROBLEMSIZE] = size(ECPDec);
    [naECP,~]          = size(ArchiveDEC);

    %% ECP optimization
    ECPDec_NEW  = [];
    ECPDec_NEW1 = [];
    ECPDec_NEW2 = [];

    for i = 1 : ceil(nECP/pop_fac)
        % generate beta = Gaussian number, with mean=0.7 and standard
        % deviation=0.2
        Force = normrnd(0.7,0.2);

        SP = sort(randperm(nECP,nECPI));

        SP = SP(1:nECPI);

        if Strategy == 1
            for ii = 1 : nECPI
                for jj = 1 : nECPI
                    S1 = ECPDec(SP(ii),:)+ Force * (ECPDec(1,:) - ECPDec(SP(ii),:));
                    if jj < ii
                        S1 = S1 + Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                        ECPDec_NEW(end+1,:) = S1;
                    elseif jj > ii
                        S1 = S1 - Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                        ECPDec_NEW(end+1,:) = S1;
                    end
                end
            end
        elseif Strategy == 2
            for ii = 1 : nECPI
                S1 = ECPDec(SP(ii),:)+ 0*Force * (ECPDec(1,:) - ECPDec(SP(ii),:));
                for jj = 1 : nECPI
                    if jj < ii
                        S1 = S1 + Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                    elseif jj > ii
                        S1 = S1 - Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                    end
                end
                ECPDec_NEW(end+1,:) = S1;
            end
        elseif Strategy == 3
            for ii = 1 : nECPI
                S2 = ECPDec(SP(ii),:)+ Force * (ECPDec(1,:) - ECPDec(SP(ii),:));
                for jj = 1 : nECPI
                    S1 = ECPDec(SP(ii),:) + Force * (ECPDec(1,:) - ECPDec(SP(ii),:));
                    if jj < ii
                        S1 = S1 + Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                        S2 = S2 + Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                        ECPDec_NEW1(end+1,:) = S1;
                    elseif jj > ii
                        S1 = S1 - Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                        S2 = S2 - Force * (ECPDec(SP(jj),:) - ECPDec(SP(ii),:));
                        ECPDec_NEW1(end+1,:) = S1;
                    end
                end
                ECPDec_NEW2(end+1,:) = S2;
            end
            ECPDec_NEW = [ECPDec_NEW1; ECPDec_NEW2];
        end
    end

    ECPDec_NEW = bound(Problem,ECPDec_NEW);

    for i1 = 1 : size(ECPDec_NEW,1)
        for j = 1 : PROBLEMSIZE
            r = rand;
            if r < 0.2
                pos = randi(naECP(1));
                ECPDec_NEW(i1,j) = ArchiveDEC(pos,j);
            end
        end
    end
    ECP_NEW = Problem.Evaluation(ECPDec_NEW);
end

function x = bound(Problem,x)
    x = min(max(x,Problem.lower),Problem.upper);
end