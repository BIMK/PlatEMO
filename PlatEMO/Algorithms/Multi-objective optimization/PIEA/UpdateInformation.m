function indicator = UpdateInformation(flag,score,indicator)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yang Li (email: liyangwust@163.com)

    a    = indicator(1).Choose_record;
    b    = indicator(2).Choose_record;
    c    = indicator(3).Choose_record;
    a(1) = [];
    b(1) = [];
    c(1) = [];
    
    if flag == 1
        a = [a 1];
        b = [b 0];
        c = [c 0];
    elseif flag == 2
        a = [a 0];
        b = [b 1];
        c = [c 0];
    else
        a = [a 0];
        b = [b 0];
        c = [c 1];
    end
    indicator(1).Choose_record = a;
    indicator(2).Choose_record = b;
    indicator(3).Choose_record = c;
    
    a = indicator(1).Win_record;
    b = indicator(2).Win_record;
    c = indicator(3).Win_record;
    
    a(1) = [];
    b(1) = [];
    c(1) = [];
    
    if score == 0
        a = [a 0];
        b = [b 0];
        c = [c 0];
    else
        if flag == 1
            a = [a score/2];
            b = [b 0];
            c = [c 0];
        elseif flag == 2
            a = [a 0];
            b = [b score/2];
            c = [c 0];
        else
            a = [a 0];
            b = [b 0];
            c = [c score/2];
        end
    end
    indicator(1).Win_record = a;
    indicator(2).Win_record = b;
    indicator(3).Win_record = c;
    
    p = [(eps+sum(indicator(1).Win_record))/(eps+sum(indicator(1).Choose_record))...
         (eps+sum(indicator(2).Win_record))/(eps+sum(indicator(2).Choose_record))...
         (eps+sum(indicator(3).Win_record))/(eps+sum(indicator(3).Choose_record))];
    
    p = p/sum(p);
    indicator(1).Pw = p(1);
    indicator(2).Pw = p(2);
    indicator(3).Pw = p(3);
end