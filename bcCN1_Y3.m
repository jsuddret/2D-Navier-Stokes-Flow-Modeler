function [a, b, c, d] = bcCN1_Y3(a, b, c, d, t)
% define global variables
global yc3;
y = yc3;
% define spacing
[M, N] = size(a); M = M-6; N = N-6;
for j = 1:N+6
    for i = 1:M+6
        % left boundaries
        if i == 4
            % inlet boundary (Dirichlet, cell-centered)
            if y(j) >= 1 && y(j) <= 1+0.75
                b(i, j) = b(i, j)-a(i, j);
                d(i, j) = d(i, j)-a(i, j)*2*Y1(t);
                a(i, j) = 0;
            % wall boundary (Neumann, cell-centered)
            else
                b(i, j) = b(i, j)+a(i, j);
                a(i, j) = 0;
            end
        end
        % right boundaries (Neumann, cell-centered)
        if i == M+3
            % outlet boundary and wall boundary 
            b(i, j) = b(i, j)+c(i, j);
            c(i, j) = 0;
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
