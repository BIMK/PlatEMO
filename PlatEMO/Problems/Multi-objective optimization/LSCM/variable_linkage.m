function [PopDec] = variable_linkage(PopDec,obj)
PPP = PopDec;
for i = 1 : obj.duan
    N = size(PopDec,1);
    if i ==  1
        P = PopDec(:,obj.CV{i,1});
        last_P = PPP(:,1);
        PP = repmat(last_P,1,size(P,2)) + P;
        a2 = mod(PP,0.5);
        if obj.CCT == 1
            PopDec(:,obj.CV{i,1}) = -2*a2+1;
        elseif obj.CCT == 2
            PopDec(:,obj.CV{i,1}) = cos(a2*pi);
        end
        control_D = obj.Dis{i,1}{1,1}(1) : obj.Dis{i,obj.M}{obj.nk(i),1}(end);
        if obj.DCT==1
            PopDec(:,control_D) = (1+repmat((1:length(control_D))./length(control_D),N,1)).*PopDec(:,control_D) - repmat(last_P,1,length(control_D));
        elseif obj.DCT==2
            PopDec(:,control_D) = (1+ cos( 0.5*pi* repmat((1:length(control_D))./length(control_D),N,1)) ).*PopDec(:,control_D) - repmat(last_P,1,length(control_D));
        end
    else
        P = PopDec(:,obj.CV{i,1});
        last_P = PPP(:, obj.Dis{i-1,obj.M}{obj.nk(i-1),1}(end) );
        PP = repmat(last_P,1,size(P,2)) + P;
        a2 = mod(PP,0.5);
        if obj.CCT == 1
            PopDec(:,obj.CV{i,1}) = -2*a2+1;
        elseif obj.CCT == 2
            PopDec(:,obj.CV{i,1}) = cos(a2*pi);
        end
        control_D = obj.Dis{i,1}{1,1}(1) : obj.Dis{i,obj.M}{obj.nk(i)}(end);
        if obj.DCT==1
            PopDec(:,control_D) = (1+repmat((1:length(control_D))./length(control_D),N,1)).*PopDec(:,control_D) - repmat(last_P,1,length(control_D));
        elseif obj.DCT==2
            PopDec(:,control_D) = (1+ cos( 0.5*pi* repmat((1:length(control_D))./length(control_D),N,1)) ).*PopDec(:,control_D) - repmat(last_P,1,length(control_D));
        end
    end
end
end

