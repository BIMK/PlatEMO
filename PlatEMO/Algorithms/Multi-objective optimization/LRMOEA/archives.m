classdef archives < handle
    
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        dec;
        mask;
        obj;
        mobj;
        mr;
        tno;
        Gn;
    end
    methods      
        function [Arch] = archives(Pop,Dec,Mask)
            if nargin > 0
                Arch(1,size(Pop,1)) = Arch;
                for i = 1 : size(Pop,1)
                    Arch(i).dec  = Dec(i,:);
                    Arch(i).mask = Mask(i,:);
                    Arch(i).obj  = Pop(i,:);
                    Arch(i).mobj = Pop(i,:);
                    Arch(i).mr   = 0;
                    Arch(i).tno  = 0;
                    Arch(i).Gn   = 1;
                end
            end
        end
        function value = decs(Arch)
            value = cat(1,Arch.dec);
        end
        function value = objs(Arch)
            value = cat(1,Arch.obj);
        end
        function value = masks(Arch)
            value = cat(1,Arch.mask);
        end
        function value = mrs(Arch)
            value = cat(1,Arch.mr);
        end
        function value = tnos(Arch)
            value = cat(1,Arch.tno);
        end
        function value =Gns(Arch)
            value = cat(1,Arch.Gn);
        end
    end
end