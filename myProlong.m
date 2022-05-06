function P = myProlong(e2h)
    % extract the size of the fine mesh error input
    [M, N] = size(e2h);
    M = M-2; N = N-2;
    eh = zeros(M*2+2, N*2+2); % preallocate with ghost cells
    for j = 2:N+1
        for i = 2:M+1
            % for 2d cell centered
            eh(2*i-2, 2*j-2) = e2h(i, j);
            eh(2*i-1, 2*j-2) = e2h(i, j);
            eh(2*i-2, 2*j-1) = e2h(i, j);
            eh(2*i-1, 2*j-1) = e2h(i, j);
        end
    end
    % apply boundary conditions and return course mesh error output
    P = bcGS(eh);
end
