function [pop,pop1,pop2,tipe] = DLCM(kneeArray1,kneeArray2,pop1,pop2,tipe)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    pietheta  = [];
    samplePop = [];
    u         = [];
    vec       = kneeArray1-kneeArray2;    
    Cw11      = mean(kneeArray1',1);  % centroid of time K-1
    Cw12      = mean(kneeArray2',1);  % centroid of time K-2
    Dw11      = Cw11'-Cw12';  % their difference    
    for ok = 1 : size(kneeArray1,2)
        bi(ok) = norm(kneeArray1(:,ok)-kneeArray2(:,ok))./norm(Dw11);
    end
    biO  = mean(bi)/size(bi,2);            
    bi   = bi+normrnd(0,biO);
    bian = abs(bi)'.*repmat(Dw11',size(kneeArray1,2),1);
    for j = 1 : size(vec,2)    
        for i = 1 : size(vec,1)-2
            kjl   = vec(i+1:end,j);
            fenzi = sqrt(sum(kjl.^2));
            pietheta(i,j) = atan(fenzi/(vec(i,j)));
        end
    end
    i = size(pietheta,1)+1;
    for j = 1 : size(vec,2)    
       kjl = vec(i+1,j);
       pietheta(i,j) = atan(kjl/(vec(i,j)));
    end
     
    for i = 1 : size(vec,1)-2
       kjl   = Dw11(i+1:end);
       fenzi = sqrt(sum(kjl.^2));
       DW_pietheta(i) = atan(fenzi/(Dw11(i)));
       if DW_pietheta(i) < 0
           DW_pietheta(i) = DW_pietheta(i)+pi;
       end
    end
    i   = size(DW_pietheta,2)+1;
    kjl = Dw11(i+1);
    DW_pietheta(i) = atan(kjl/(Dw11(i)));     
    
    vecb = vec';
    for num = 1 : size(kneeArray1,2)+1
        if num ~= size(kneeArray1,2)+1
            fi = reshape(pietheta(:,num),[size(pietheta,1) 1]);
            for i = 1 : size(vec,1)
                if i == 1
                    u(i) = norm(vecb(num,:))*cos(fi(i));
                elseif i < size(vec,1)
                    temp = 1;
                    for j = 1 : i-1
                        temp = temp*sin(fi(j));
                    end
                    u(i) = norm(vecb(num,:))*temp*cos(fi(i));
                else
                    temp = 1;
                    for j = 1 : i-1
                        temp = temp*sin(fi(j));
                    end
                    u(i) = norm(vecb(num,:))*temp;
                end
            end
            samplePop(:,num) = u;
        else
            fi = reshape(DW_pietheta,[size(pietheta,1) 1]);
            for i = 1 : size(vec,1)
                if i == 1
                    u(i) = norm(Dw11)*cos(fi(i));
                elseif i < size(vec,1)
                    temp = 1;
                    for j = 1 : i-1
                        temp = temp*sin(fi(j));
                    end
                    u(i) = norm(Dw11)*temp*cos(fi(i));
                else
                    temp = 1;
                    for j = 1 : i-1
                        temp = temp*sin(fi(j));
                    end
                    u(i) = norm(Dw11)*temp;
                end
            end
            D412 = u;
        end      
    end    
    pam    = sum(kneeArray2+samplePop-kneeArray1,1);
    D41pam = sum(Cw12+D412-Cw11); 
    for k = 1 : size(pam,2)
        if abs(pam(k)) > 1
            if pietheta(end,k) > 0
                pietheta(end,k) = pietheta(end,k)-pi;
            else
                pietheta(end,k) = pietheta(end,k)+pi;
            end
        end
    end
    
    if abs(D41pam) > 1
        if DW_pietheta(end) > 0
            DW_pietheta(end) = DW_pietheta(end)-pi;
        else
            DW_pietheta(end) = DW_pietheta(end)+pi;
        end
    end
      
    samplePop = [];
    u = [];
    for num = 1 : size(kneeArray1,2)
        fi = reshape(pietheta(:,num),[size(pietheta,1) 1]);
        fi = 1/2*(DW_pietheta'+fi);
        for i = 1 : size(vec,1)
            if i == 1
                u(i) = norm(bian(num,:))*cos(fi(i));
            elseif i < size(vec,1)
                temp = 1;
                for j = 1 : i-1
                    temp = temp*sin(fi(j));
                end
                u(i) = norm(bian(num,:))*temp*cos(fi(i));
            else
                temp = 1;
                for j = 1 : i-1
                    temp = temp*sin(fi(j));
                end
                u(i) = norm(bian(num,:))*temp;
            end
        end
        samplePop(:,num) = u;
    end    
    pop = mod(abs(kneeArray1+samplePop),1)'; %DCM
    if tipe == 0
        selecta = size(pop,1)/2;
        tipe    = selecta;
    else
        selecta = tipe;
    end 
    CP_gl   = selecta/size(pop,1);
    CP_num  = find(rand(1,size(pop,1))<CP_gl);
    DCM_num = setdiff(1:size(pop,1),CP_num);
    pop(CP_num,:) = kneeArray1(:,CP_num)'+bian(CP_num,:); %CP
    pop1 = [pop1;pop(CP_num,:)];
    pop2 = [pop2;pop(DCM_num,:)];
end