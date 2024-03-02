function [out1,out2,CCT,DCT,Dis_function,Q_form,PF_form] = duan_yueshu(index)

if index ==1
    nk = [5     4     6     6     3     5     5     5     5     2]; % the distance function variables in each set is divided into nk subsets
    c_index = [2, 2, 2, 2 ,2, 2, 2 ,2, 2, 2];  %  CEC2006 constraint index used in this function
    CCT = 1; % Constraint variable linkage: CCT has two kinds values (1 : linear, 2 : non-linear)
    DCT = 2; % Unconstraint variable linkage: DCT has two kinds values (1 : linear, 2 : non-linear)
    Dis_function = [1,2]; % distance function index used in this function
    PF_form = 1; % PF_form has three kinds of forms (1 : linear，2 : Convex，3 : disconnected)
    Q_form = 1; % Q_form has three kinds of forms (1 : single，2 : multiple，3 : full)
    
elseif index == 2
    nk = [5     4     6     6     3     5     5     5     5     2];
    c_index = [4, 4, 4, 4 ,4, 4, 4 ,4, 4, 4]; 
    CCT = 1;
    DCT = 1;
    Dis_function = [1,1];
    PF_form = 3;
    Q_form = 1;
    
elseif index == 3
    nk = [5     4     6     6     3     5     5     5     5     2];
    c_index = [10, 10, 10, 10 ,10, 10, 10 ,10, 10, 10];   
    CCT = 1;
    DCT = 1;
    Dis_function = [1,3];
    PF_form = 2; 
    Q_form = 2;
    
elseif index == 4
    nk = [5     4    3    5      4    3    5     4     3    5];
    c_index = [4, 11, 4, 11, 4, 11, 4, 11, 4, 11]; 
    CCT = 2;
    DCT = 2;
    Dis_function = [4,3];
    PF_form = 3; 
    Q_form = 1;
    
elseif index == 5
    nk = [5     4    3    5    4    3   5    4    3    5];
    c_index = [2, 1, 2, 1, 2, 1, 2, 1, 2, 1]; 
    CCT = 1;
    DCT = 1;
    Dis_function = [2,3];
    PF_form = 1; 
    Q_form = 3;
    
elseif index == 6
    nk = [5     3     5     4     3     5     4    3     5     2];
    c_index = [7, 8, 7, 8, 7, 8, 7, 8, 7, 8]; 
    CCT = 2;
    DCT = 1;
    Dis_function = [1,5];
    PF_form = 1; 
    Q_form = 1;
    
elseif index == 7  
    nk = [3     6     5     6     3     6     5     6     3     6];
    c_index = [3, 11, 4, 3, 11, 4, 3, 11, 4, 3]; 
    CCT = 1;
    DCT = 1;
    Dis_function = [4,5];
    PF_form = 1; 
    Q_form = 2;
    
elseif index == 8
    nk = [4     4     5    6     3     5     5     5     5     2];
    c_index = [10, 11, 9, 10, 11, 9, 10, 11, 9, 10]; 
    CCT = 2;
    DCT = 2;
    Dis_function = [1,1];
    PF_form = 2; 
    Q_form = 1;
    
elseif index == 9
    nk = [4, 5, 4, 5, 4, 5, 4, 5, 4, 5];
    c_index = [4, 7, 11, 4, 7, 11, 4, 7, 11, 4]; 
    CCT = 1;
    DCT = 1;
    Dis_function = [1,4];
    PF_form = 2; 
    Q_form = 2;
    
elseif index == 10
    nk = [4, 5, 4, 5, 4, 5, 4, 5, 4, 5];
    c_index = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3];
    CCT = 2;
    DCT = 2;
    Dis_function = [1,2];
    PF_form = 1; 
    Q_form = 1;
    
elseif index == 11
    c_index = [4, 6, 4, 6, 4, 6, 4, 6, 4, 6, 4, 6];
    nk = [2, 3, 2, 3, 2, 3, 2, 3, 2, 3];
    CCT = 1;
    DCT = 2;
    Dis_function = [2,3];
    PF_form = 2; 
    Q_form = 2;
    
elseif index == 12
    c_index = [4,5,6,4,4,5,6,4,4,5];
    nk = [11, 4, 3, 11, 4, 3, 11, 4, 3, 11];
    CCT = 1;
    DCT = 1;
    Dis_function = [1,1];
    PF_form = 3; 
    Q_form = 1;
end

CEC_problem = [1,2,3,6,9,10,11,12,14,18,24]; %  CEC2006 constraint problem index used in this study

out1 = nk;
out2 = CEC_problem(c_index);
end

