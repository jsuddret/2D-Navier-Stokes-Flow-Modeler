function P = myPoisson(phi, f, h, nIterMax, epsilon)
    %{
    * : including ghost cells
    phi : 2d array of cell centered variable with initial guess *
    f : 2d array of cellc entered PDE RHS values *
    h : equidistant mesh spacing in the x and y-direction
    nIterMax : integer number of maximum V-cycle iterations
    epsilon : convergence threshold for infinity norm of the residual
    %}
    phi = bcGS(phi);
    iterator = 0;
    error = 10^6;
    while iterator < nIterMax && error >= epsilon
        phi = myMultigrid(phi, f, h);
        res = myResidual(phi, f, h); % computer residual of phi
        res = reshape(res, [], 1); % convert matrix to row vector
        error = max(abs(res));
        iterator = iterator+1;
    end
    % return 2d array of cell centered variable with solution *
    P = bcGS(phi);
end