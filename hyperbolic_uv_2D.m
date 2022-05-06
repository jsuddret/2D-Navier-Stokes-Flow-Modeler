function [Hu, Hv] = hyperbolic_uv_2D(u, v)
% define global variables
global h; % step size
% determine M, N from u
[M, N] = size(u); M = M-1; N = N-2; % size of interior
% preallocate outputs
Hu = zeros(size(u));
Hv = zeros(size(v));
% calculate partials for u
u1 = V(u,2:M,2:N+1,1:M-1,2:N+1);
u2 = V(u,3:M+1,2:N+1,2:M,2:N+1);
duudx = (u2.^2-u1.^2)/h;
u1 = V(u,2:M,2:N+1,2:M,3:N+2).*V(v,2:M,2:N+1,3:M+1,2:N+1);
u2 = V(u,2:M,1:N,2:M,2:N+1).*V(v,2:M,1:N,3:M+1,1:N);
duvdy = (u1-u2)/h;
% populate Hu
Hu(2:M,2:N+1) = -(duudx+duvdy);
% calcualte partials for v
v1 = V(v,2:M+1,2:N,2:M+1,1:N-1);
v2 = V(v,2:M+1,3:N+1,2:M+1,2:N);
duudx = (v2.^2-v1.^2)/h;
v1 = V(u,2:M+1,2:N,2:M+1,3:N+1).*V(v,2:M+1,2:N,3:M+2,2:N);
v2 = V(u,1:M,2:N,1:M,3:N+1).*V(v,1:M,2:N,2:M+1,2:N);
duvdy = (v1-v2)/h;
% populate Hv
Hv(2:M+1,2:N) = -(duudx+duvdy);
end

% 2nd-order central
function velocity = V(vel,i1,j1,i2,j2)
velocity = 1/2*(vel(i1,j1)+vel(i2,j2));
end