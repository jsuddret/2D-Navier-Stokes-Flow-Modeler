function p = parabolic_CN1_2D_v(v, Qu, dt)
% define global variables
global a b c d;
global Re t h;
% define sizing
[M, N] = size(v); M = M-2; N = N-1;
% define coefficients
dx = 1/Re*dt/h^2;
d1 = dx/2;
d2 = d1;
% define a, b, c, d vectors with the same size as u
m = ones(M+2, N+1); % temporary value to create vectors
a = -d1.*m; b = (1+2*d1).*m; c = a; d = m;
% loop through interior cells and calculate d + source term at t^n
for j = 2:N % cell-centered interior
    for i = 2:M+1 % node-based interior
        d(i, j) = d2*v(i, j+1)+(1-2*d2)*v(i, j)+d2*v(i, j-1)+Qu(i, j)*dt/2;
    end
end
% apply implicit boundary conditions
[a, b, c, d] = bcCN1_v(a, b, c, d, t+dt/2);
% reshape from 2D to 1D
a = reshape(a(2:M+1, 2:N), [(N-1)*M, 1]);
b = reshape(b(2:M+1, 2:N), [(N-1)*M, 1]);
c = reshape(c(2:M+1, 2:N), [(N-1)*M, 1]);
d = reshape(d(2:M+1, 2:N), [(N-1)*M, 1]);
% solve tridiagonal matrix and reshape to 2D
v(2:M+1, 2:N) = reshape(mySolveTriDiag(a, b, c, d), [M, N-1]);
% reshape from 1D to 2D and store a, b, c, d
tempa = reshape(a, [M, N-1]); tempb = reshape(b, [M, N-1]);
tempc = reshape(c, [M, N-1]); tempd = reshape(d, [M, N-1]);
% create as the same size as u
a = zeros(M+2, N+1); a(2:M+1, 2:N) = tempa;
b = zeros(M+2, N+1); b(2:M+1, 2:N) = tempb;
c = zeros(M+2, N+1); c(2:M+1, 2:N) = tempc;
d = zeros(M+2, N+1); d(2:M+1, 2:N) = tempd;
% apply boundary conditions to u^(n+1/2)
v = bc_v(v, t+dt/2);
% return u^(n+1/2)
p = v;
end