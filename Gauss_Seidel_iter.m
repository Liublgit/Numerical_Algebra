function x = Gauss_Seidel_iter(A,b,tol,it_max)
% Gauss_Seidel_iter solve the problem Ax = b
% using Gauss-Seidel iteration method.
% tol : norm(x_k-x_{k-1}) < tol
% it_max : max of iteration.
%=============== initiate  ===============%
iter = 1;
n = size(A);
x0 = ones(n(1),1);
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
B = (D - L)\U;
g = (D - L)\b;
rho = vrho(B);
fprintf('\n Spectral radius of iteration matrix = %d \n',rho)
%=============== iteration  ===============%
while iter <= it_max
    x = B*x0 + g;
    fprintf('\n Iteration number = %d, ||x_k - x_{k-1}||= %d \n',iter, norm(x-x0))
    if norm(x-x0) < tol
        break;
    end
    x0 = x;
    iter = iter + 1;
end
%=============== Exceed  ===============%
if iter > it_max
   fprintf('\n Maximum number of iterations exceeded! \n');
end
end
