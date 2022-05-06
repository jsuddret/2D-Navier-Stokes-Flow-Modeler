function w = calc_adYdx_WENO_2D(Y,a)
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
dp(1:M+5,4:N+3) = Y(2:M+6,4:N+3)-Y(1:M+5,4:N+3);
dmdp(2:M+5,4:N+3) = Y(3:M+6,4:N+3)-2*Y(2:M+5,4:N+3)+Y(1:M+4,4:N+3);
% loop through interior
for j = 4:N+3
    for i = 4:M+3
        if a(i,j) <= 0
            WENOa = dmdp(i+2,j)/h;
            WENOb = dmdp(i+1,j)/h;
            WENOc = dmdp(i,j)/h;
            WENOd = dmdp(i-1,j)/h;
            f(i,j) = a(i,j)*(1/(12*h)*(-dp(i-2,j)+7*dp(i-1,j)+...
                7*dp(i,j)-dp(i+1,j))+...
                psiWENO(WENOa, WENOb, WENOc, WENOd));
        elseif a(i,j) > 0
            WENOa = dmdp(i-2,j)/h;
            WENOb = dmdp(i-1,j)/h;
            WENOc = dmdp(i,j)/h;
            WENOd = dmdp(i+1,j)/h;
            f(i,j) = a(i,j)*(1/(12*h)*(-dp(i-2,j)+7*dp(i-1,j)+...
                7*dp(i,j)-dp(i+1,j))-...
                psiWENO(WENOa, WENOb, WENOc, WENOd));
        end
    end
end
% return f
w = f;
end