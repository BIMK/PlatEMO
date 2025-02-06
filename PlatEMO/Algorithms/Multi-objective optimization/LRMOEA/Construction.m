function [rf,mf,ar,wr] = Construction(arch,W,score)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [~,n]  = size(arch);
    [N ,~] = size(W);
    if size(arch,2) > size(W,1)
        a(N+1).r = [];
        de = arch.decs;
        ma = arch.masks;
        m  = arch.Gns;
        A_sp = zeros(1,n);
        Sp   = zeros(1,n);
        for i = 1 : n
            a(arch(i).tno).r = [a(arch(i).tno).r,i];
            sp   = sum(arch(i).mobj,2);
            sp_b = min(sp);
            A_sp(1,i) = mean(sum(abs(sp - sp_b)));
            A_sp(1,i) = 1/((sum(ma(i,:)&(score>0.5)))/size(score,2));
            Sp(1,i)   = mean(sp);
        end
        rch = [];
        w   = [];
        j   = 0;
        for r = 1 : N
            if a(r).r > 0
                [~,c] = size(a(r).r);    
                s     = zeros(1,c);
                t     = zeros(1,c);
                gn    = 0;
                for d = 1 : c
                    gn = gn + a(r).r(1,d);
                end
                gn = gn/c;
                for d = 1 : c
                    no = a(r).r(1,d);
                    if( m(no)>gn )
                        T = 0;
                    else
                        T = gn-m(no);
                    end
                    s(1,d) = Sp(1, no)+(A_sp(1,no));
                    t(1,d) = no;
                end
                [~,b] = sort(s(1,:));
                j  = j + 1;
                tt = t(1,b(1,1));
                RF(j,:) = de(t(1,b(1,1)),:);
                MF(j,:) = ma(t(1,b(1,1)),:);
                rch = [rch,tt];
                w   = [w,r];
            end
        end
    end
    ar = rch;
    wr = w;
    rf = RF;
    mf = MF;
end