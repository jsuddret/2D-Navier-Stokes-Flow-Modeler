function H = hyperbolic_Y_WENO_2D(Y,u0,v0,u1,v1,dt)
% define global variables
global t; % current time unit
global a0 b0 a1 b1 a2 b2;
global Y1 Y2 Ys;
% determine mesh size from Y
[M, N] = size(Y);
% remove ghost cells
M = M-6; N = N-6;
% cell-center u velocities
a0 = zeros(size(Y)); a1 = zeros(size(Y));
a0(4:M+3,4:N+3) = 1/2*(u0(2:M+1,2:N+1)+u0(1:M,2:N+1));
a1(4:M+3,4:N+3) = 1/2*(u1(2:M+1,2:N+1)+u1(1:M,2:N+1));
% second order u velocity at t^(n+1/2)
a2 = (a1+a0)/2;
% cell-center v velocities
b0 = zeros(size(Y)); b1 = zeros(size(Y));
b0(4:M+3,4:N+3) = 1/2*(v0(2:M+1,2:N+1)+v0(2:M+1,1:N));
b1(4:M+3,4:N+3) = 1/2*(v1(2:M+1,2:N+1)+v1(2:M+1,1:N));
% second order v velocity at t^(n+1/2)
b2 = (b1+b0)/2;
% step 1
Y1 = Y-dt*(calc_adYdx_WENO_2D(Y,a0)+calc_bdYdy_WENO_2D(Y,b0));
Y1 = bc_Y3(Y1,t+dt);
% step 2
Y2 = Y1+3/4*dt*(calc_adYdx_WENO_2D(Y,a0)+calc_bdYdy_WENO_2D(Y,b0))-...
    1/4*dt*(calc_adYdx_WENO_2D(Y1,a1)+calc_bdYdy_WENO_2D(Y1,b1));
Y2 = bc_Y3(Y2,t+1/2*dt);
% step 3
Y3 = Y2+1/12*dt*(calc_adYdx_WENO_2D(Y,a0)+calc_bdYdy_WENO_2D(Y,b0))+...
    1/12*dt*(calc_adYdx_WENO_2D(Y1,a1)+calc_bdYdy_WENO_2D(Y1,b1))-...
    2/3*dt*(calc_adYdx_WENO_2D(Y2,a2)+calc_bdYdy_WENO_2D(Y2,b2));
Ys = bc_Y3(Y3,t+dt);
HY = (Ys-Y)/dt;
H = HY;
end