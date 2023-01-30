function P = MatingSelection(S)
% Mating selection of C-TAEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
  
    number=length(S);
    rnd=randi(number,1,2);
    x_1=rnd(1);
    x_2=rnd(2);
    CV1 = sum(max(0,S(x_1).con),2);
    CV2 = sum(max(0,S(x_2).con),2);
    if CV1>0 && CV2>0
        x=randi(number,1);
        P=S(x);
    elseif CV1<=0 && CV2>0
        P=S(x_1);
    elseif CV1>0  && CV2<=0
        P=S(x_2);
    elseif CV1==0 && CV2==0
          [FrontNo,~] = NDSort(S(rnd).objs,inf);
          if FrontNo(1)<=FrontNo(2)
              P=S(x_1);
          else
              P=S(x_2);
          end
     end
end