function orthonormal_basis = Gram_Schmidt_process(A)
%
% Gram-Schmidt Process
% 
% Gram_Schmidt_process(A) produces an orthonormal basis for the subspace of
% Eucldiean n-space spanned by the vectors {u1,u2,...}, where the matrix A 
% is formed from these vectors as its columns. That is, the subspace is the
% column space of A. The columns of the matrix that is returned are the 
% orthonormal basis vectors for the subspace. An error is returned if an
% orthonormal basis for the zero vector space is attempted to be produced.
%
% For example, if the vector space V = span{u1,u2,...}, where u1,u2,... are
% row vectors, then set A to be [u1' u2' ...].
%
% For example, if the vector space V = Col(B), where B is an m x n matrix,
% then set A to be equal to B.

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

        for i = 2:n
            u = A(1:m,i);
            v = zeros(m,1);

            for j = 1:(i - 1)
                v = v - dot(u,orthonormal_basis(1:m,j))*orthonormal_basis(1:m,j);
            end

            v_ = u + v;
            orthonormal_basis(1:m,i) = v_/norm(v_);
        end
    end
end