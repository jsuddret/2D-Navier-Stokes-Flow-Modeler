function p = parabolic_CN2_2D_Y3(Y, QY, dt)
% define global variables
global a b c d;
global Re Sc t h;
% define sizing
[M, N] = size(Y); M = M-6; N = N-6;
% define coefficients
dx = 1/(Re*Sc)*dt/h^2;
d1 = dx/2;
d2 = d1;
% define a, b, c, d vectors with the same size as u
m = ones(M+6, N+6); % temporary value to create vectors
a = -d2.*m; b = (1+2*d2).*m; c = a; d = m;
% loop through interior cells and calculate d + source term at t^n
for j = 4:N+3 % cell-centered interior
    for i = 4:M+3 % node-based interior
        d(i, j) = d1*Y(i+1, j)+(1-2*d1)*Y(i, j)+d1*Y(i-1, j)+QY(i, j)*dt/2;
    end
end
% apply implicit boundary conditions
[a, b, c, d] = bcCN2_Y3(a, b, c, d, t+dt);
% reshape from 2D to 1D
a = reshape(a(4:M+3, 4:N+3)', [N*M, 1]);
b = reshape(b(4:M+3, 4:N+3)', [N*M, 1]);
c = reshape(c(4:M+3, 4:N+3)', [N*M, 1]);
d = reshape(d(4:M+3, 4:N+3)', [N*M, 1]);
% solve tridiagonal matrix and reshape to 2D
Y(4:M+3, 4:N+3) = reshape(mySolveTriDiag(a, b, c, d), [N, M])';
% reshape from 1D to 2D and store a, b, c, d
tempa = reshape(a, [N, M])'; tempb = reshape(b, [N, M])';
tempc = reshape(c, [N, M])'; tempd = reshape(d, [N, M])';
% create as the same size as u
a = zeros(M+6, N+6); a(4:M+3, 4:N+3) = tempa;
b = zeros(M+6, N+6); b(4:M+3, 4:N+3) = tempb;
c = zeros(M+6, N+6); c(4:M+3, 4:N+3) = tempc;
d = zeros(M+6, N+6); d(4:M+3, 4:N+3) = tempd;
% apply boundary conditions to u^(n+1/2)
Y = bc_Y3(Y, t+dt);
% return u^(n+1/2)
p = Y;
end