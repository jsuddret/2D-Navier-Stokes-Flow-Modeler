function bc = bcGhost_u(u, t)
% define global variables
global xf;
% define local variable from global
x = xf;
% define spacing
[M, N] = size(u); M = M-1; N = N-2;
for j = 1:N+2
    for i = 1:M+1
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
% return u
bc = u;
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