function result = is_orthogonal_set(A)
%
% Orthogonal Sets
% 
% is_orthogonal_set(A) determines if a set of vectors in Euclidean n-space 
% is orthogonal. The matrix A is formed from these vectors as its columns. 
% That is, the subspace spanned by the set of vectors is the column space 
% of A. The value 1 is returned if the set is orthogonal. The value 0 is 
% returned if the set is not orthogonal.
%
% For example, if the set of row vectors (u1,u2,...} is to be determined 
% for orthogonality, set A to be equal to [u1' u2' ...].

matrix_size = size(A);

m = matrix_size(1,1);
n = matrix_size(1,2);

tolerance = 10^-10;

if n == 1
        result = 1;
else
    orthogonal_counter = 0;

    for i = 1:n
        for j = 1:n
            if i == j
            else
                if abs(dot(A(1:m,i),A(1:m,j))) <= tolerance
                    orthogonal_counter = orthogonal_counter + 1;
                end
            end
        end
    end

    if orthogonal_counter == factorial(n)/factorial(n - 2);
        result = 1;
    else
        result = 0;
    end
end