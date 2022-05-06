function M = myMultigrid(phi, f, h)
    %{
    * : including ghost cells
    phi : 2d array of a cell centered solution variable *
    f : 2d arrray of the cell centered PDE RHS *
    h : equidistant mesh space in the x and y-direction
    %}
    [M, N] = size(phi); % extract the size of phi
    M = M-2; N = N-2; % eliminate ghost cells from M, N values
    n = 1; % perform a single V-cycle
    nc = 4; % four additional iterations on the coursest mesh
    phi = myGaussSeidel(phi, f, h, n);
    % base case 
    if mod(M, 2)==0 && mod(N, 2)==0 % M, N must be divisible by two
        rh = myResidual(phi, f, h);
        r2h = myRestrict(rh); % restrict mesh
        e2h = zeros(M/2+2, N/2+2); % preallocate 2d cell centered mesh
        e2h = myMultigrid(e2h, r2h, 2*h); % recursion call
        eh = myProlong(e2h); % prolong mesh
        phi = phi+eh;
        phi = bcGS(phi); % apply boundary conditions
        phi = myGaussSeidel(phi, f, h, n); % apply gauss-seidel
    else
       phi = myGaussSeidel(phi, f, h, nc); % apply gauss-seidel
    end
    % apply boundary conditions and return course mesh error output
    M = bcGS(phi);
end