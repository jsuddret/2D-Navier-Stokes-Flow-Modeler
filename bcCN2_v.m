function [a, b, c, d] = bcCN2_v(a, b, c, d, t)
% define global variables
global xc;
x = xc;
% define spacing
[M, N] = size(a); M = M-2; N = N-1;
for j = 1:N+1
    for i = 1:M+2
        % bottom boundaries (Dirichlet)
        if j == 2
            % inlet boundary
            if x(i) >= 1.5 && x(i) <= 1.5+0.75
                d(i, j) = d(i, j)-a(i, j)*v2(x(i), t);
                a(i, j) = 0; % neighbors do not exist in linear system
            % wall boundary
            else
                a(i, j) = 0; % neighbors do not exist in linear system
            end
        end
        % top boundaries (Dirichlet)
        if j == N
            if x(i) >= 0.75 && x(i) <= 0.75+0.5
            % inlet boundary
                d(i, j) = d(i, j)-c(i, j)*v3(x(i));
                c(i, j) = 0; % neighbors do not exist in linear system
            % wall boundary
            else
                c(i, j) = 0; % neighbors do not exist in linear system
            end
        end
    end
end
end

% v2 function
function v = v2(x, t)
z = z2(x);
f = fi(z);
v = f*sin(a2(t));
    function a = a2(t)
        a = pi/2+pi/4*sin(2*pi/5*t);
    end
    function z = z2(x)
        z = (x-1.5)/0.75;
    end
end

% v3 function
function v = v3(x)
z = z3(x);
f = fi(z);
v = -f;
    function z = z3(x)
        z = (x-0.75)/0.5;
    end
end

function f = fi(z)
f = 6*z*(1-z);
end
