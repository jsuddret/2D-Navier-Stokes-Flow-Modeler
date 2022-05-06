function d = calcDivV(u, v)
% define global variables
global h;
% determine M, N from u
[M, N] = size(u); M = M-1; N = N-2;
% calculate partial u
du = zeros(M+2, N+2);
du(2:M+1, 2:N+1) = (u(2:M+1, 2:N+1)-u(1:M, 2:N+1))/h;
% calculate partial v
dv = zeros(M+2, N+2);
dv(2:M+1, 2:N+1) = (v(2:M+1, 2:N+1)-v(2:M+1, 1:N))/h;
% create divV
divV = du+dv;
% return divV
d = divV;
end