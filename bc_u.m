function bc = bc_u(u, t)
% define global variables
global xf yc;
x = xf; y = yc;
% define spacing
[M, N] = size(u); M = M-1; N = N-2;
% apply boundary conditions
for j = 1:N+2
    for i = 1:M+1
        % left boundaries (Dirichlet, node-based)
        if i == 1
            % inlet boundary
            if y(j) >= 1 && y(j) <= 1+0.75
                u(i, j) = u1(y(j), t);
            % wall boundary
            else
                u(i, j) = 0;
            end
        end
        % right boundaries
        if i == M+1
            % outlet boundary (Neumann, node-based)
            if y(j) >= 0.25 && y(j) <= 0.25+1.25
                u(i, j) = (-1*u(i-2, j)+4*u(i-1, j))/3;
            % wall boundary (Dirichlet, node-based)
            else
                u(i, j) = 0;
            end
        end
        % bottom boundaries (Dirichlet, cell-centered)
        if j == 1
            % inlet boundary
            if x(i) >= 1.5 && x(i) <= 1.5+0.75
                u(i, j) = 2*u2(x(i), t)-u(i, j+1);
            % wall boundary
            else
                u(i, j) = -u(i, j+1);
            end
        end
        % top boundaries (Dirichlet, cell-centered)
        if j == N+2
            % inlet boundary and wall boundary
            u(i, j) = -u(i, j-1);
        end
    end
end
% return u with boundary conditions applied
bc = u;
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
