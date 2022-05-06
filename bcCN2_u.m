function [a, b, c, d] = bcCN2_u(a, b, c, d, t)
% define global variables
global xf;
x = xf;
% define spacing
[M, N] = size(a); M = M-1; N = N-2;
for j = 1:N+2
    for i = 1:M+1
        % bottom boundaries (Dirichlet)
        if j == 2
            % inlet boundary
            if x(i) >= 1.5 && x(i) <= 1.5+0.75
                b(i, j) = b(i, j)-a(i, j);
                d(i, j) = d(i, j)-a(i, j)*2*u2(x(i), t);
                a(i, j) = 0;
            % wall boundary
            else
                b(i, j) = b(i, j)-a(i, j);
                a(i, j) = 0;
            end
        end
        % top boundaries (Dirichlet)
        if j == N+1
            % inlet boundary and wall boundary
            b(i, j) = b(i, j)-c(i, j);
            c(i, j) = 0;
        end
    end
end
end

% u2 function
function u = u2(x, t)
z = z2(x);
f = fi(z);
u = f*cos(a2(t));
    function a = a2(t)
        a = pi/2+pi/4*sin(2*pi/5*t);
    end
    function z = z2(x)
        z = (x-1.5)/0.75;
    end
end

function f = fi(z)
f = 6*z*(1-z);
end
