function FS = Recombination(AC,AD,Xic,Xre,eps_peak,fS,RankC)
% Recombination

%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    ACDec = AC.decs;
    ADDec = AD.decs;     
    [N,D] = size(ADDec);
    peak  = find( ((fS - min(fS)) < eps_peak) & RankC == 1 );    
    FS    = NaN(length(peak).*N,D);
    FS(:,Xic) =  kron( ACDec(peak',Xic),ones(N,1));
    FS(:,Xre) =  repmat(ADDec(:,Xre),length(peak),1);          
end