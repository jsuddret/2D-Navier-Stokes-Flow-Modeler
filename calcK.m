function c = calcK(u, v)
% define global variables
global Lx Ly h;
% attain M and N from size of u
[M, N] = size(u);
% remove ghost cells
M = M-1; N = N-2;
% average u velocity
uavg = zeros(M, N+2);
for j = 1:N+2
    for i = 1:M
        uavg(i, j) = 1/2*(u(i, j)+u(i+1, j));
    end
end
% remove ghost cells
uavg = uavg(1:M, 2:N+1);
% avgerage v velocity
vavg = zeros(M+2, N);
for i = 1:M+2
    for j = 1:N
        vavg(i, j) = 1/2*(v(i, j)+v(i, j+1));
    end
end
% remove ghost cells
vavg = vavg(2:M+1, 1:N);
% compute average velocity over interior
avg = 1/2.*(uavg.^2+vavg.^2);
% initialize k
k = 0;
% double summation
for j = 1:N
    for i = 1:M
        k = k+avg(i, j)*h^2;
    end
end
% return S
c = 1/(Lx*Ly)*k;
end