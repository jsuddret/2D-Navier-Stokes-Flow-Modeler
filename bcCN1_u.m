function [a, b, c, d] = bcCN1_u(a, b, c, d, t)
% define global variables
global yc;
y = yc;
% define spacing
[M, N] = size(a); M = M-1; N = N-2;
for j = 1:N+1
    for i = 1:M+1
        % left boundaries (Dirichlet)
        if i == 2
            % inlet boundary
            if y(j) >= 1 && y(j) <= 1+0.75
                d(i, j) = d(i, j)-a(i, j)*u1(y(j), t);
                a(i, j) = 0; % neighbors do not exist in linear system
            % wall boundary
            else
                a(i, j) = 0; % neighbors do not exist in linear system
            end
        end
        % right boundaries
        if i == M
            % outlet boundary (Neumann)
            if y(j) >= 0.25 && y(j) <= 0.25+1.25
                b(i, j) = b(i, j)+4/3*c(i, j);
                a(i, j) = a(i, j)-1/3*c(i, j);
                c(i, j) = 0; % neighbors do not exist in linear system
            % wall boundary (Dirichlet)
            else
                c(i, j) = 0; % neighbors do not exist in linear system
            end
        end
    end
end
end

% u1 function
function u = u1(y, t)
z = z1(y);
f = fi(z);
u = f*cos(a1(t));
    function a = a1(t)
        a = pi/3*sin(pi*t);
    end
    function z = z1(y)
        z = (y-1)/0.75;
    end
end

function f = fi(z)
f = 6*z*(1-z);
end
