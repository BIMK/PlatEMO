function [winner,loser] = CHP(P1,P2,epsilon) 
% Paired competition

% Copyright (c) 2020-2021 Cheng He

    winner = P1;
    loser  = P2;
    N      = length(P1);
    CV1    = sum(max(P1.cons,0),2);
    CV2    = sum(max(P2.cons,0),2);
    flag   = zeros(N,1);
    for i = 1 : N
        if max(CV1(i),CV2(i)) <= epsilon
            flag(i) = 2 - (norm(P1(i).objs)<norm(P2(i).objs));
        else
            flag(i) = 2 - (CV1(i)<CV2(i));
        end
        if flag(i) == 1
			[winner(i),loser(i)] = deal(P1(i),P2(i));
        else
            [winner(i),loser(i)] = deal(P2(i),P1(i));
        end
    end
end