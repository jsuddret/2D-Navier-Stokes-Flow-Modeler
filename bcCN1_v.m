function [a, b, c, d] = bcCN1_v(a, b, c, d, t)
% define global variables
global yf;
y = yf;
% define spacing
[M, N] = size(a); M = M-1; N = N-2;
for j = 1:N+2
    for i = 1:M+1
        % left boundaries (Dirichlet, cell-centered)
        if i == 2
            % inlet boundary
            if y(j) >= 1 && y(j) <= 1+0.75
                b(i, j) = b(i, j)-a(i, j);
                d(i, j) = d(i, j)-a(i, j)*2*v1(y(j), t);
                a(i, j) = 0;
            % wall boundary
            else
                b(i, j) = b(i, j)-a(i, j);
                a(i, j) = 0;
            end
        end
        % right boundaries
        if i == M
            % outlet boundary (Neumann, cell-centered)
            if y(j) >= 0.25 && y(j) <= 0.25+1.25
                b(i, j) = b(i, j)+c(i, j);
                c(i, j) = 0;
            % wall boundary (Dirichlet, cell-centered)
            else
                b(i, j) = b(i, j)-c(i, j);
                c(i, j) = 0;
            end
        end
    end
end
end

% v1 function
function v = v1(y, t)
z = z1(y);
f = fi(z);
v = f*sin(a1(t));
    function a = a1(t)
        a = pi/3*sin(pi*t);
    end
    function z = z1(y)
        z = (y-1)/0.75;
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