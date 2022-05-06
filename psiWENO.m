function p = psiWENO(a,b,c,d)
% define epsilon
epsilon = 10^-6;
% define ISi terms
IS0 = 13*(a-b).^2+3*(a-3*b).^2;
IS1 = 13*(b-c).^2+3*(b+c).^2;
IS2 = 13*(c-d).^2+3*(3*c-d).^2;
% define alpha terms
alpha0 = 1./(epsilon+IS0).^2;
alpha1 = 6./(epsilon+IS1).^2;
alpha2 = 3./(epsilon+IS2).^2;
% define omega terms
omega0 = alpha0./(alpha0+alpha1+alpha2);
omega2 = alpha2./(alpha0+alpha1+alpha2);
% calculate psi
psi = 1/3*omega0.*(a-2*b+c)+1/6*(omega2-1/2).*(b-2*c+d);
% return psi
p = psi;
end