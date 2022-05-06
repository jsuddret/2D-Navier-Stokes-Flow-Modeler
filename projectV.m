function [u, v] = projectV(u, v, phi, dt)
% define global variables
global t h;
% find mesh size from phi
[M, N] = size(phi); M = M-2; N = N-2;
% preallocate the gradient of phi
dphidx = zeros(size(phi));
dphidy = zeros(size(phi));
% calculate the gradient of phi
dphidx(2:M+1, 2:N+1) = (phi(3:M+2, 2:N+1)-phi(2:M+1, 2:N+1))/h;
dphidy(2:M+1, 2:N+1) = (phi(2:M+1, 3:N+2)-phi(2:M+1, 2:N+1))/h;
% project onto u, v and apply boundary conditions at next time step
u(2:M, 2:N+1) = u(2:M, 2:N+1)-dt*dphidx(2:M, 2:N+1);
v(2:M+1, 2:N) = v(2:M+1, 2:N)-dt*dphidy(2:M+1, 2:N);
% apply ghost cell boundary conditions to u, v
u = bcGhost_u(u, t+dt);
v = bcGhost_v(v, t+dt);
end
