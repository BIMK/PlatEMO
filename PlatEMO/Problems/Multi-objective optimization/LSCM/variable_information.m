function [obj] = variable_information(obj,c)
% Dis{i,j}{h}，i indicates the index of obj.duan，j indicates the index of objectives，h indicates the index of nk
for j = 1 : obj.duan
    if j == 1
        obj.sublen{j} = floor(c./sum(c).*(obj.base - obj.M - obj.high_D_C(j))/obj.nk(j));
        obj.len{j}    = [0,cumsum(obj.sublen{j}*obj.nk(j))];
        obj.CV{j,1} =  obj.M +1 :  obj.M + obj.high_D_C(j);
        
        temp = obj.len{j} + obj.CV{j}(end);
        for M = 1 : obj.M
            for i = 1 : obj.nk(j)
                obj.Dis{j,M}{i,1} = temp(M)+1 + (i-1)* obj.sublen{j}(M):  temp(M) + i*obj.sublen{j}(M);
            end
        end
    else 
        obj.sublen{j} = floor(c./sum(c).*(obj.base - obj.high_D_C(j))/obj.nk(j));
        obj.len{j}    = [0,cumsum(obj.sublen{j}*obj.nk(j))];
        obj.CV{j,1} =  (j-1)*obj.base + 1 : (j-1)*obj.base + obj.high_D_C(j);
        
        temp = obj.len{j} + obj.CV{j}(end);
        for M = 1 : obj.M
            for i = 1 : obj.nk(j)
                obj.Dis{j,M}{i,1} = temp(M)+1 + (i-1)* obj.sublen{j}(M):  temp(M) + i*obj.sublen{j}(M);
            end
        end
    end
end
end

