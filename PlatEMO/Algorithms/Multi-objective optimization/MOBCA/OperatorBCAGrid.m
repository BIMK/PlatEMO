function Offspring = OperatorBCAGrid(Problem,Parent,BCB,Archive,nSoldiers,nArmies)
%The operator of MOBCA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    ParticleDec = Parent.decs;
    ArchiveDec  = Archive.decs;
    [N,D]       = size(ParticleDec);
    Offspring   = [];
    
    %% Particle swarm optimization
    for i=1:nArmies
        for j=1:nSoldiers
            r=randi(nArmies,1);
            while(i==r)
                r=randi(nArmies,1);
            end
            for d=1:D
                if rand<BCB
                    alpha=rand*2*pi;
                    soldiers(j,d)=ArchiveDec(i,d)+abs(ParticleDec(r,d)-ParticleDec(i,d))*sin(alpha);
                else
                    beta=rand*2*pi;
                    soldiers(j,d)=ParticleDec(r,d)+abs(ParticleDec(r,d)-ParticleDec(i,d))*cos(beta);
                end% end BCB
                soldiers(j,d)=max(min(soldiers(j,d),Problem.upper(d)),Problem.lower(d));
            end% end D
        end% end nSoldiers
        Offspring=[Offspring;soldiers];
    end% end nArmies
    Offspring = Problem.Evaluation(Offspring);
end