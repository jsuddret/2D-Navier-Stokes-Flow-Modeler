function [u, v] = correctOutlet(u, v)
% define global variables
global yc h ucorr;
% set global variable to local variable
y = yc;
% determine M, N from u
[M, N] = size(u); M = M-1; N = N-2;
% initialize Lout based on schematic
Lout = 1.25;
% initialize qdot* to zero
qdotstar = 0;
% loop through j in u interior
for j = 2:N+1
    qdotstar = qdotstar+u(1, j)*h-u(M+1, j)*h;
end
% loop through i in v interior
for i = 2:M+1
    qdotstar = qdotstar+v(i, 1)*h-v(i, N+1)*h;
end
% determine correction velocity
ucorr = qdotstar/Lout;
% add ucorr to right outlet
% loop through u interior
for j = 2:N+1
    if y(j) >= 0.25 && y(j) <= 0.25+1.25
        u(M+1, j) = u(M+1, j)+ucorr;
    end
end
% return u, v
end