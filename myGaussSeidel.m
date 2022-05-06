function GS = myGaussSeidel(phi, f, h, niter)
%{
    * : including ghost cells
    phi : 2d array of cell centered PDE solution variable *
    f : 2d array of cell centered PDE RHS *
    h : mesh spacing
    niter : integer number of iterations
[M, N] = size(phi); % number of points in [x, y] *
for k = 1:niter % loop through number of iterations
    phi(2:M-1, 2:N-1) = 0.25*(phi(1:M-2, 2:N-1)+phi(3:M, 2:N-1)+...
        phi(2:M-1, 1:N-2)+phi(2:M-1, 3:N))-0.25*h^2*f(2:M-1, 2:N-1);
    phi = bcGS(phi); % apply boundary conditions
end
GS = phi; % return solution
%}
[M, N] = size(phi); % number of points in [x, y] *
for k = 1:niter % loop through number of iterations
    for j = 2:N-1 % loop through stencil points in y
        for i = 2:M-1 % loop through stencil points in x
            phi(i, j) = 0.25*(phi(i-1, j)+phi(i+1, j)+...
                phi(i, j-1)+phi(i, j+1))-0.25*h^2*f(i, j);
        end
    end
    phi = bcGS(phi); % apply boundary conditions
end
GS = phi; % return solution
end
