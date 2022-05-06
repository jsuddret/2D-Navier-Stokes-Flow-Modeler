function R = myRestrict(rh)
    [M, N] = size(rh); % extract the size of the fine mesh residual input
    M = (M-2)/2; N = (N-2)/2; % divide interior mesh by two
    r2h = zeros(M+2, N+2); % preallocate with ghost cells
    for j = 2:N+1 % loop through stencil points in y
       for i = 2:M+1 % loop through stencil points in x
          r2h(i, j) = 1/4*(rh(2*i-2, 2*j-2)+rh(2*i-1, 2*j-2)+...
              rh(2*i-2, 2*j-1)+rh(2*i-1, 2*j-1)); % for 2d cell centered
       end
    end
    R = r2h; % return course mesh residual output
end
