function bc = bc_Y3(Y, t)
% define global variables
global xc3 yc3;
x = xc3; y = yc3;
% define spacing
[M, N] = size(Y); M = M-6; N = N-6;
% apply boundary conditions
for j = 1:N+6
    for i = 1:M+6
        % left boundaries 
        if i == 3
            % inlet boundary (Dirichlet, cell-centered)
            if y(j) >= 1 && y(j) <= 1+0.75
                % apply to all three ghost cells
                for temp = 1:3
                    Y(i-(temp-1), j) = 2*Y1(t)-Y(i+temp, j);
                end
            % wall boundary (Neumann, cell-centered)
            else
                % apply to all three ghost cells
                for temp = 1:3
                   Y(i-(temp-1), j) =  Y(i+temp, j);
                end
            end
        end
        % right boundaries
        if i == M+4
            % outlet boundary and wall boundary (Neumann, cell-centered)
            for temp = 1:3
                Y(i+(temp-1), j) = Y(i-temp, j);
            end
        end
        % bottom boundaries
        if j == 3
            % inlet boundary (Dirichlet, cell-centered)
            if x(i) >= 1.5 && x(i) <= 1.5+0.75
                % apply to all three ghost cells
                for temp = 1:3
                    Y(i, j-(temp-1)) = 2*Y2(t)-Y(i, j+temp);
                end
            % wall boundary (Neumann, cell-centered)
            else
                % apply to all three ghost cells
                for temp = 1:3
                   Y(i, j-(temp-1)) =  Y(i, j+temp);
                end
            end
        end
        % top boundaries
        if j == N+4
            if x(i) >= 0.75 && x(i) <= 0.75+0.5
            % inlet boundary (Dirichlet, cell-centered)
                % apply to all three ghost cells
                for temp = 1:3
                   Y(i, j+(temp-1)) = 2*Y3(t)-Y(i, j-temp);
                end
            % wall boundary (Neumann, cell-centered)
            else
                % apply to all three ghost cells
                for temp = 1:3
                   Y(i, j+(temp-1)) = Y(i, j-temp);
                end
            end
        end
    end
end
% return Y with boundary conditions applied
bc = Y;
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
