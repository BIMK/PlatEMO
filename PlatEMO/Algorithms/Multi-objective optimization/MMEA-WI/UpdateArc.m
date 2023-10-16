function [Population,dk] = UpdateArc(Population,offspring,N)
% Update the archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenhua Li

    joint=[Population offspring];

    [FrontNo,MaxFNo] = NDSort(joint.objs,N);
    next=FrontNo==1;
    Population=joint(next);

    [Choose,fDKN] = DoubleNearestSelection(Population,N,5);
    Population = Population(Choose);
    dk = fDKN(Choose);
end

function [Choose,fDN] = DoubleNearestSelection(Pop,N,K)
    PopObj=Pop.objs;
    PopDec=Pop.decs;
    Np = size(PopObj,1);
    Choose = true(1,Np);
    if Np <= K
        fDN = zeros(1,Np);
        return;
    end
    d_obj = pdist2(PopObj,PopObj,'euclidean');
    d_dec = pdist2(PopDec,PopDec,'euclidean');
    d_obj(logical(eye(Np))) = inf;
    d_dec(logical(eye(Np))) = inf;

    sdo = sort(d_obj);
    sdd = sort(d_dec);
    dn_obj = sum(sdo(1:K,:));
    dn_dec = sum(sdd(1:K,:));
    avg_dn_obj = mean(dn_obj);
    avg_dn_dec = mean(dn_dec);
    if avg_dn_obj == 0
        avg_dn_obj = inf;
    end
    if avg_dn_dec == 0
        avg_dn_dec = inf;
    end
    fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);

    while sum(Choose) > N
        [~,Del] = max(fDN);

        Choose(Del) = false;
        d_obj(Del,:) = inf;
        d_obj(:,Del) = inf;
        d_obj(Del,:) = inf;
        d_obj(:,Del) = inf;

        sdo = sort(d_obj);
        sdd = sort(d_dec);
        dn_obj = sum(sdo(1:K,:));
        dn_dec = sum(sdd(1:K,:));
        avg_dn_obj = mean(dn_obj);
        avg_dn_dec = mean(dn_dec);
        if avg_dn_obj == 0
            avg_dn_obj = inf;
        end
        if avg_dn_dec == 0
            avg_dn_dec = inf;
        end
        fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
        fDN(~Choose) = -inf;
    end
end