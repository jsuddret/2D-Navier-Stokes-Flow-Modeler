function [a, b, c, d] = bcCN2_Y3(a, b, c, d, t)
% define global variables
global xc3;
x = xc3;
% define spacing
[M, N] = size(a); M = M-6; N = N-6;
for j = 1:N+6
    for i = 1:M+6
        % bottom boundaries
        if j == 4
            % inlet boundary (Dirichlet, cell-centered)
            if x(i) >= 1.5 && x(i) <= 1.5+0.75
                b(i, j) = b(i, j)-a(i, j);
                d(i, j) = d(i, j)-a(i, j)*2*Y2(t);
                a(i, j) = 0;
            % wall boundary (Neumann, cell-centered)
            else
                b(i, j) = b(i, j)+a(i, j);
                a(i, j) = 0;
            end
        end
        % top boundaries
        if j == N+3
            if x(i) >= 0.75 && x(i) <= 0.75+0.5
            % inlet boundary (Dirichlet, cell-centered)
                b(i, j) = b(i, j)-c(i, j);
                d(i, j) = d(i, j)-c(i, j)*2*Y3(t);
                c(i, j) = 0;
            % wall boundary (Neumann, cell-centered)
            else
                b(i, j) = b(i, j)+c(i, j);
                c(i, j) = 0;
            end
        end
    end
end
end

function y = Y1(t)
y = fY(t);
end

function y = Y2(t)
y = 1-fY(t);
end

function y = Y3(t)
y = 1-fY(t);
end

function f = fY(t)
if mod(t, 4) < 2
    f = 1;
else
    f = 0;
end
end
