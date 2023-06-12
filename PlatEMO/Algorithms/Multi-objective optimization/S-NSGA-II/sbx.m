function Offspring = sbx(Parent, lb, ub, Parameter)
% Migrated from the PlatEMO OperatorGA module for convenience

% This function is written by Ian Meyer Kropp

    [proC,disC] = deal(Parameter{:});

    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);

    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    
    Lower = repmat(lb,2*N,1);
    Upper = repmat(ub,2*N,1);
    
    % Put everything back in bounds
    Offspring = min(max(Offspring,Lower),Upper);
end