function [k,d,mu,lm] = qpsubp(dfk,Bk,Ae,hk,Ai,gk)
% Solve the quadratic programming subproblem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    n = length(dfk);
    l = length(hk);
    m = length(gk);
    beta    = 0.5;
    sigma   = 0.2;
    epsilon = 1e-6;
    gamma   = 0.05;
    ep0     = 0.05;
    d0      = ones(n,1);
    mu0     = 0.05*zeros(l,1);
    lm0     = 0.05*zeros(m,1);
    u0      = [ep0;zeros(n+l+m,1)];
    z0 = [ep0;d0;mu0;lm0];
    z  = z0;
    ep = ep0;
    d  = d0;
    mu = mu0;
    lm = lm0;
    for k = 1 : 150
        dh = dah(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk);
        if norm(dh) < epsilon
            break;
        end
        A  = JacobiH(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk);
        b  = psi(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk,gamma)*u0 - dh;
        dz = A\b;
        de = dz(1);
        dd = dz(2:n+1);
        if l>0 && m>0
            du = dz(n+2:n+l+1);
            dl = dz(n+l+2:n+l+m+1);
        elseif l == 0
            dl = dz(n+2:n+m+1);
        elseif m == 0
            du = dz(n+2:n+l+1);
        end
        for mk = 0 : 20
            t1 = beta^mk;
            if l>0 && m>0
                dh1 = dah(ep+t1*de,d+t1*dd,mu+t1*du,lm+t1*dl,dfk,Bk,Ae,hk,Ai,gk);
            elseif l == 0
                dh1 = dah(ep+t1*de,d+t1*dd,mu,lm+t1*dl,dfk,Bk,Ae,hk,Ai,gk);
            elseif m == 0
                dh1 = dah(ep+t1*de,d+t1*dd,mu+t1*du,lm,dfk,Bk,Ae,hk,Ai,gk);
            end
            if norm(dh1) <= (1-sigma*(1-gamma*ep0)*beta^mk)*norm(dh)
                break;
            end
        end
        alpha = beta^mk;
        ep    = ep + alpha*de;
        d     = d + alpha*dd;
        if l>0 && m>0
            mu = mu + alpha*du;
            lm = lm + alpha*dl;
        elseif l == 0
            lm = lm + alpha*dl;
        elseif m == 0
            mu = mu + alpha*du;
        end
    end
end

function p = phi(ep,a,b)
    p = a + b - sqrt(a^2+b^2+2*ep^2);
end

function dh = dah(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk)
    n  = length(dfk);
    l  = length(hk);
    m  = length(gk);
    dh = zeros(n+l+m+1,1);
    dh(1) = ep;
    if l>0 && m>0
        dh(2:n+1) = Bk*d - Ae'*mu - Ai'*lm + dfk;
        dh(n+2:n+l+1) = hk + Ae*d;
        for i = 1 : m
            dh(n+l+1+i) = phi(ep,lm(i),gk(i)+Ai(i,:)*d);
        end
    elseif l == 0
        dh(2:n+1) = Bk*d - Ai'*lm + dfk;
        for i = 1 : m
            dh(n+1+i) = phi(ep,lm(i),gk(i)+Ai(i,:)*d);
        end
    elseif m == 0
        dh(2:n+1) = Bk*d - Ae'*mu + dfk;
        dh(n+2:n+l+1) = hk + Ae*d;
    end
    dh = dh(:);
end

function xi = psi(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk,gamma)
    dh = dah(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk);
    xi = gamma*norm(dh)*min(1,norm(dh));
end

function [dd1,dd2,v1] = ddv(ep,d,lm,Ai,gk)
    m   = length(gk);
    dd1 = zeros(m);
    dd2 = zeros(m);
    v1  = zeros(m,1);
    for i = 1 : m
        fm = sqrt(lm(i)^2+(gk(i)+Ai(i,:)*d)^2+2*ep^2);
        dd1(i,i) = 1 - lm(i)/fm;
        dd2(i,i) = 1 - (gk(i)+Ai(i,:)*d)/fm;
        v1(i) = -2*ep/fm;
    end
end

function A = JacobiH(ep,d,mu,lm,dfk,Bk,Ae,hk,Ai,gk)
    n = length(dfk);
    l = length(hk);
    m = length(gk);
    [dd1,dd2,v1] = ddv(ep,d,lm,Ai,gk);
    if l>0 && m>0
        A = [1          zeros(1,n) zeros(1,l) zeros(1,m)
             zeros(n,1) Bk         -Ae'       -Ai'
             zeros(l,1) Ae         zeros(l,l) zeros(l,m)
             v1         dd2*Ai     zeros(m,l) dd1];
    elseif l == 0
        A = [1          zeros(1,n) zeros(1,m)
             zeros(n,1) Bk         -Ai'
             v1         dd2*Ai     dd1];
    elseif m == 0
        A = [1          zeros(1,n) zeros(1,l)
             zeros(n,1) Bk         -Ae'
             zeros(l,1) Ae         zeros(l)];
    end
end