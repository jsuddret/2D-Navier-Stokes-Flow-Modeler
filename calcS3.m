function c = calcS3(Y)
% define global variables
global Lx Ly h;
% attain M and N from size of input
[M, N] = size(Y);
% remove ghost cells
M = M-6; N = N-6;
% initialize S
S = 0;
% double summation
for j = 4:N+3
    for i = 4:M+3
        S = S+R(Y(i, j))*h^2; % dA=h*h=h^2
    end
end
% return S
c = 1/(Lx*Ly).*S;
end

% function call for R as a function of Y
function r = R(Y)
r = Y*(1-Y);
end