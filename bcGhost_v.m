function bc = bcGhost_v(v, t)
% define global variables
global yf;
% define local variable from global
y = yf;
% define spacing
[M, N] = size(v); M = M-2; N = N-1;
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
    end
    % return v
    bc = v;
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

function f = fi(z)
f = 6*z*(1-z);
end