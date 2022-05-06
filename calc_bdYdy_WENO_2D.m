function w = calc_bdYdy_WENO_2D(Y,a)
% define global variables
global h; % step size
% determine mesh size from Y
[M, N] = size(Y);
% remove ghost cells
M = M-6; N = N-6;
% define f from size Y
f = zeros(size(Y)); % set ghost layers to zero
% define dp, dpms
dp = zeros(size(Y));
dmdp = zeros(size(Y));
% populate dp, dpms
dp(4:M+3,1:N+5) = Y(4:M+3,2:N+6)-Y(4:M+3,1:N+5);
dmdp(4:M+3,2:N+5) = Y(4:M+3,3:N+6)-2*Y(4:M+3,2:N+5)+Y(4:M+3,1:N+4);
% loop through interior
for j = 4:N+3
    for i = 4:M+3
        if a(i,j) <= 0
            WENOa = dmdp(i,j+2)/h;
            WENOb = dmdp(i,j+1)/h;
            WENOc = dmdp(i,j)/h;
            WENOd = dmdp(i,j-1)/h;
            f(i,j) = a(i,j)*(1/(12*h)*(-dp(i,j-2)+7*dp(i,j-1)+...
                7*dp(i,j)-dp(i,j+1))+...
                psiWENO(WENOa, WENOb, WENOc, WENOd));
        elseif a(i,j) > 0
            WENOa = dmdp(i,j-2)/h;
            WENOb = dmdp(i,j-1)/h;
            WENOc = dmdp(i,j)/h;
            WENOd = dmdp(i,j+1)/h;
            f(i,j) = a(i,j)*(1/(12*h)*(-dp(i,j-2)+7*dp(i,j-1)+...
                7*dp(i,j)-dp(i,j+1))-...
                psiWENO(WENOa, WENOb, WENOc, WENOd));
        end
    end
end
% return f
w = f;
end