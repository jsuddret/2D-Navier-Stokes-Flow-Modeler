% solves a tridiagonal linear system using Gaussian elimination
function sol = mySolveTriDiag(a, b, c, d)
    % error handling
    try
        % if input argumments are not of identical lengths, 
        % vertcat will throw an exception
        temp = [a; b; c; d];
        % dummy function to handle warning
        unused(temp);
    catch
        error("Input arguments must be of identical length.");
    end
    
    % elimination
    P = length(a); % number of rows in matrix
    for i = 2:P % loop through rows, exclude first row
        b(i) = b(i) - c(i-1)*a(i)/b(i-1);
        d(i) = d(i) - d(i-1)*a(i)/b(i-1);
    end
    
    % back substitution
    d(P) = d(P)/b(P);
    for i = P-1:-1:1 % loop through rows backwards, exclusing last row
       d(i) = (d(i) - c(i)*d(i+1))/b(i);
    end
    
    % return solution
    sol = d;
end

% prevents unused warning at line 5
function t = unused(temp)
    t = temp;
end
