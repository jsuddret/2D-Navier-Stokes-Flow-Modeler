function [dt, outputFlag] = calcDt561(t, outputTime, u, v)
% define global variables
global h CFL;
% define a, b, c, d
a = 2*u;
b = v;
c = u;
d = 2*v;
if sum(sum(a)) == 0 && sum(sum(b)) == 0 || sum(sum(c)) == 0 && sum(sum(d)) == 0
    dt = 1e-6;
else
    % define time step for u and v
    dthu = CFL*h/(max(max(abs(a)))+max(max(abs(b))));
    dthv = CFL*h/(max(max(abs(c)))+max(max(abs(d))));
    % find smaller of the two time steps
    dt = min([dthu, dthv]);
end
% initialize output flag to zero
outputFlag = 0;
if t < outputTime && t + dt >= outputTime
    dt = outputTime-t;
    % since dt was adjusted, set to 1
    outputFlag = 1;
end
end