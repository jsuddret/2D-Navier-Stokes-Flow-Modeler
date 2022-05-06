function r = myResidual(phi, f, h)
%{
    * : including ghost cells
    phi : 2d array of cell centered PDE solution variable *
    f : 2d array of cell centered PDE RHS *
    h : mesh spacing
%}
[M, N] = size(phi); % number of points in [x, y] *
res = zeros(M, N); % preallocate residual vector
res(2:M-1, 2:N-1) = f(2:M-1, 2:N-1)-...
    1/h^2*((phi(1:M-2, 2:N-1)-2*phi(2:M-1, 2:N-1)+phi(3:M, 2:N-1))+...
    (phi(2:M-1, 1:N-2)-2*phi(2:M-1, 2:N-1)+phi(2:M-1, 3:N)));
r = res; % return residual vector
end