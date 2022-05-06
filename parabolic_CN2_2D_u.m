function p = parabolic_CN2_2D_u(u, Qu, dt)
% define global variables
global a b c d;
global Re t h;
% define sizing
[M, N] = size(u); M = M-1; N = N-2;
% define coefficients
dx = 1/Re*dt/h^2;
d1 = dx/2;
d2 = d1;
% define a, b, c, d vectors with the same size as u
m = ones(M+1, N+2); % temporary value to create vectors
a = -d2.*m; b = (1+2*d2).*m; c = a; d = m;
% loop through interior cells and calculate d + source term at t^n
for j = 2:N+1 % cell-centered interior
    for i = 2:M % node-based interior
        d(i, j) = d1*u(i+1, j)+(1-2*d1)*u(i, j)+d1*u(i-1, j)+Qu(i, j)*dt/2;
    end
end
% apply implicit boundary conditions
[a, b, c, d] = bcCN2_u(a, b, c, d, t+dt);
% reshape to turn 2D arrays into 1D
a = reshape(a(2:M, 2:N+1)', [(M-1)*N, 1]);
b = reshape(b(2:M, 2:N+1)', [(M-1)*N, 1]);
c = reshape(c(2:M, 2:N+1)', [(M-1)*N, 1]);
d = reshape(d(2:M, 2:N+1)', [(M-1)*N, 1]);
% solve tridiagonal matrix and reshape to 2D
u(2:M, 2:N+1) = reshape(mySolveTriDiag(a, b, c, d), [N, M-1])';
% reshape from 1D to 2D and store a, b, c, d
tempa = reshape(a, [N, M-1])'; tempb = reshape(b, [N, M-1])';
tempc = reshape(c, [N, M-1])'; tempd = reshape(d, [N, M-1])';
% create as the same size as u
a = zeros(M+1, N+2); a(2:M, 2:N+1) = tempa;
b = zeros(M+1, N+2); b(2:M, 2:N+1) = tempb;
c = zeros(M+1, N+2); c(2:M, 2:N+1) = tempc;
d = zeros(M+1, N+2); d(2:M, 2:N+1) = tempd;
% apply boundary conditions to u^(n+1/2)
u = bc_u(u, t+dt);
% return u^(n+1/2)
p = u;
end