function OffDec = ACO_MV_generate(Parent,len_r,len_c,k,l,rank_v,m,up_r,dn_r)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    q    = 0.05099;
    kesi = 0.6795;
    w    = ( 1/(q*k*sqrt(2*pi)) )*exp( -((rank_v-1).^2)/(2*(q*k)^2) );
    p_ro = w/sum(w);
    v_r  = Parent(:,1:len_r);
    v_c  = Parent(:,len_r+1:len_r+len_c);
    for i = 1 : m
        v_r_generate(i,:) = ACO_RO(v_r,p_ro,len_r,k,kesi);
        v_r_generate(i,:) = Repair(v_r_generate(i,:),up_r,dn_r);
        v_c_generate(i,:) = ACO_C(v_c,l,len_c,q,w);
    end
    OffDec = [v_r_generate,v_c_generate];
end
    
function [n_x] = ACO_C(x,l,len_c,q,w)
    for j = 1 : len_c
        pl       = Cal_pl(x(:,j),l(:,j),w,q);
        idx_gc   = Select(pl);
        n_x(1,j) = idx_gc;
    end
end
    
function [out] = ACO_RO(x,p,len,k,kesi)
    global zzz BBB v000 xxxxx idx_grrr
    idx_gr   = Select(p);
    [z,B,v0] = Rot(x,idx_gr,len);
    idx_grrr = idx_gr;
    zzz      = z;
    BBB      = B;
    v000     = v0;
    xxxxx    = x;
    for j = 1 : len
        mu       = z(idx_gr,j);
        sigma    = kesi*sum( abs( z(idx_gr,j)-z(:,j) ) )/(k-1);
        n_z(1,j) = mu + sigma*randn(1);
    end
    out = n_z*B'+v0;
end
    
function [z_r,B,v0] = Rot(v_r,idx_gr,len_r)
    if (sum(sum(v_r - v_r(idx_gr,:))) ~= 0) && (len_r>1)
        B = VCH(v_r,v_r(idx_gr,:));
    else
        B = eye(len_r);
    end
    if rank(B) ~= len_r
        B = eye(len_r);
    end
    z_r = (v_r - v_r(idx_gr,:))*B;
    v0  = v_r(idx_gr,:);
end
    
function i = Select(p)
    p_sel = cumsum(p);
    R     = rand;
    i     = 1;
    while p_sel(i) < R
        i = i + 1;
    end
end
    
function out = Repair(x,up,dn)
    x   = (x>=up).*max(dn,2*up-x) + (x<up).*x;
    x   = (x<=dn).*min(up,2*dn-x) + (x>dn).*x;
    out = x;
end
    
function [out] = Cal_pl(vc,l,w,q)
    for i = 1 : l
        idx_l = ( vc==i );
        u(i)  = sum( idx_l );
        if isempty(w(idx_l))
            wjl(i) = 0;
        else
            wjl(i) = max(w(idx_l));
        end
    end
    eta = 100*sum( u==0 );
    for i = 1 : l
        if (eta>0)&&(u(i)>0)
            wl(i) = wjl(i)/u(i) + q/eta;
        elseif (eta==0)&&(u(i)>0)
            wl(i) = wjl(i)/u(i);
        elseif (eta>0)&&(u(i)==0)
            wl(i) = q/eta;
        end
    end
    out = wl/sum(wl);
end
    
function out = VCH(s,sl)
    [~,n] = size(s);
    for i = 1 : n
        ds       = sqrt( sum( (sl(:,i:n) - s(:,i:n)).^2, 2) );
        p        = ds.^4 / sum(ds.^4);
        idx      = Select(p);
        A(i,:)   = sl - s(idx,:);
        s(idx,:) = [];
    end
    if max(max(A)) < 1e-5
        B = Gram_Schmidt_process(rand(n));
    else
        B = Gram_Schmidt_process(A');
    end
    out = B;
end
    
function orthonormal_basis = Gram_Schmidt_process(A)
    matrix_size = size(A);
    m = matrix_size(1,1);
    n = matrix_size(1,2);
    if A == zeros(m,n)
        error('There does not exist any type of basis for the zero vector space.');
    elseif n == 1
        orthonormal_basis = A(1:m,1)/norm(A(1:m,1));
    else
        flag = 0;
        if is_orthonormal_set(A) == 1
            orthonormal_basis = A;
            flag = 1;
        end
        if flag == 0
            if rank(A) ~= n
                A = basis_col(A);
            end
            matrix_size = size(A);
            m = matrix_size(1,1);
            n = matrix_size(1,2);
            orthonormal_basis = A(1:m,1)/norm(A(1:m,1));
            for i = 2 : n
                u = A(1:m,i);
                v = zeros(m,1);
                for j = 1 : (i-1)
                    v = v - dot(u,orthonormal_basis(1:m,j))*orthonormal_basis(1:m,j);
                end
                v_ = u + v;
                orthonormal_basis(1:m,i) = v_/norm(v_);
            end
        end
    end
end
    
function basis = basis_col(A)
    matrix_size = size(A);
    m = matrix_size(1,1);
    n = matrix_size(1,2);
    if A == zeros(m,n)
        error('There does not exist a basis for the zero vector space.');
    elseif n == 1
        basis = A;
    else
        flag = 0;
        if n == 2
            multiple = A(1,2)/A(1,1);
            count    = 0;
            for i = 1 : m
                if A(i,2)/A(i,1) == multiple
                    count = count + 1;
                end
            end
            if count == m
                basis = A(1:m,1);
                flag = 1;
            end
        end
        if flag == 0
            [~,pivot_columns] = ref(A);
            for i = 1 : size(pivot_columns,2)
                B(1:m,i) = A(1:m,pivot_columns(1,i));
            end
            basis = B;
        end
    end
end
    
function [U,pivcol,nonpivcol] = ref(A,tol)
    [m,n] = size(A);
    tiny  = max(m,n)*eps*max(1,norm(A,'inf'))*10;   % default tolerance for zeros
    if (nargin==2)
        tiny = tol;                                 % reset tolerance, if specified
    end
    pivcol    = [];
    nonpivcol = [];
    U = A;
    i = 1;                                          % row index
    j = 1;                                          % column index
    while (i <= m) && (j <= n)
        [x,k] = max(abs(U(i:m,j))); p = k+i-1;      % value and row index of next pivot.
        if (x <= tiny)                              % This column is negligible.
            U(i:m,j)  = zeros(m-i+1,1);             % So clean up the entries.
            nonpivcol = [nonpivcol j];
            j = j + 1;                              % Pivot row p must be recalculated.
        else
            U([i p],j:n) = U([p i],j:n);            % Swap the ith and pth rows.
            U(i,j:n)     = U(i,j:n)/U(i,j);         % Divide pivot row by the pivot.
            for k = [1:i-1 i+1:m]                   % Replacement operations on other rows.
                U(k,j:n) = U(k,j:n) - U(k,j)*U(i,j:n);
            end
            pivcol = [pivcol j];
            i = i + 1;
            j = j + 1;
        end
    end
    nonpivcol = [nonpivcol j:n];
end
    
function result = is_orthonormal_set(A)
    matrix_size = size(A);
    m = matrix_size(1,1);
    n = matrix_size(1,2);
    tolerance = 10^-10;
    if A == zeros(m,n)
        error('The set that contains just zero vectors cannot be orthonormal.');
    elseif n == 1
        if norm(A(1:m,1)) == 1
            result = 1;
        else
            result = 0;
        end
    else
        if is_orthogonal_set(A) == 1
            length_counter = 0;
            for i = 1 : n
                if abs(norm(A(1:m,i)) - 1) <= tolerance
                    length_counter = length_counter + 1;
                end
            end
            if length_counter == n
                result = 1;
            else
                result = 0;
            end
        else
            result = 0;
        end
    end
end
    
function result = is_orthogonal_set(A)
    matrix_size = size(A);
    m = matrix_size(1,1);
    n = matrix_size(1,2);
    tolerance = 10^-10;
    if n == 1
        result = 1;
    else
        orthogonal_counter = 0;
        for i = 1 : n
            for j = 1 : n
                if i == j
                else
                    if abs(dot(A(1:m,i),A(1:m,j))) <= tolerance
                        orthogonal_counter = orthogonal_counter + 1;
                    end
                end
            end
        end
        if orthogonal_counter == factorial(n)/factorial(n - 2)
            result = 1;
        else
            result = 0;
        end
    end
end