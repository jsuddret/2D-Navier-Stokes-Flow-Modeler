function bc = bcGS(phi)
%{
* : including ghost cells
phi : 2d array of cell centered PDE solution variable *
%}
global Morig;
Mf = Morig;
[M,N] = size(phi); % number of points in [x, y] *
% left boundary conditions
phi(2-1,2:N-1) = phi(2,2:N-1);
% right boundary conditions
phi(M-1+1,2:N-1) = phi(M-1,2:N-1);
% top boundary conditions
phi(2:M-1,2-1) = phi(2:M-1,2);
% bottom boundary conditions
if M == Mf+2
    phi(2:M-1,N-1+1) = 4-phi(2:M-1,N-1);
else
    phi(2:M-1,N-1+1) = -phi(2:M-1,N-1);
end
bc = phi; % return applied boundary conditions
end