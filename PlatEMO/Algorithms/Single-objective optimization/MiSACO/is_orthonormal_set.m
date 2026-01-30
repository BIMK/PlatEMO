function result = is_orthonormal_set(A)
%
% Orthonormal Sets
% 
% is_orthonormal_set(A) determines if a set of vectors in Euclidean n-space
% is orthonormal. The matrix A is formed from these vectors as its columns. 
% That is, the subspace spanned by the set of vectors is the column space 
% of A. The value 1 is returned if the set is orthonormal. The value 0 is 
% returned if the set is not orthonormal. An error is returned if a set
% containing only zero vectors is attempted to be determined for
% orthonormality.
%
% For example, if the set of row vectors (u1,u2,...} is to be determined 
% for orthonormality, set A to be equal to [u1' u2' ...].

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
        
        for i = 1:n
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