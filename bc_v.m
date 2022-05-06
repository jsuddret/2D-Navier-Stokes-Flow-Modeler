function bc = bc_v(v, t)
% define global variables
global xc yf;
x = xc; y = yf;
% define spacing
[M, N] = size(v); M = M-2; N = N-1;
% apply boundary conditions
for j = 1:N+1
    for i = 1:M+2
        % left boundaries (Dirichlet, cell-centered)
        if i == 1
            % inlet boundary
            if y(j) >= 1 && y(j) <= 1+0.75
                v(i, j) = 2*v1(y(j), t)-v(i+1, j);
            % wall boundary
            else
                v(i, j) = -v(i+1, j);
            end
        end
        % right boundaries
        if i == M+2
            % outlet boundary (Neumann, cell-centered)
            if y(j) >= 0.25 && y(j) <= 0.25+1.25
                v(i, j) = v(i-1, j);
            % wall boundary (Dirichlet, cell-centered)
            else
                v(i, j) = -v(i-1, j);
            end
        end
        % bottom boundaries (Dirichlet, node-based)
        if j == 1
            % inlet boundary
            if x(i) >= 1.5 && x(i) <= 1.5+0.75
                v(i, j) = v2(x(i), t);
            % wall boundary
            else
                v(i, j) = 0;
            end
        end
        % top boundaries (Dirichlet, node-based)
        if j == N+1
            if x(i) >= 0.75 && x(i) <= 0.75+0.5
            % inlet boundary
                v(i, j) = v3(x(i));
            % wall boundary
            else
                v(i, j) = 0;
            end
        end
    end
end
% return u with boundary conditions applied
bc = v;
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
