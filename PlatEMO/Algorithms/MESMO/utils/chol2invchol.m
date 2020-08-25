function MatrixInv = chol2invchol(Matrix)
u = chol(Matrix);
tu=invutri(u);
MatrixInv = tu*tu';
end

function ut=invutri(u)
   [m,n]=size(u);
   if m~=n
       error('Error using in invutri.Matrix must be square.');
   end
   for i=1:n
       for j=1:i-1
           if u(i,j)~=0
               error('Error using in invutri.Matrix must be up-triangle.');
           end
       end
   end
   for i=1:n
       if u(i,i)==0
           error('Error using in invutri.Matrix is singular.')
       end
   end
   mu=[u,eye(n)];
   for i=n:-1:1
       mu(i,:)=mu(i,:)/mu(i,i);
       for j=i-1:-1:1
           mu(j,:)=mu(j,:)-mu(i,:)*mu(j,i);
       end
   end
   ut=mu(:,n+1:2*n);
end