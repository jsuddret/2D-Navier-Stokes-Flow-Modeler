clc; clear all;

% Final Project
% Jake Suddreth
% 1214553743

global CFL h Lx Ly t Re Sc xc xc3 xf yc yc3 yf Morig;

M = 64; % number of cells in the x
Morig = M;
N = 32; % number of cells in the y

% create M,N vector for GCI
M =  [M];
N = [N];
g = length(M);

% prealloce T,KE
T = zeros(g,1);
KE = zeros(g,1);

% preallocate movies
movieindex = 1;
umovie = VideoWriter('u.avi');
umovie.FrameRate = 10;
vmovie = VideoWriter('v.avi');
vmovie.FrameRate = 10;
Ymovie = VideoWriter('Y.avi');
Ymovie.FrameRate = 10;
Rmovie = VideoWriter('R.avi');
Rmovie.FrameRate = 10;

for i = 1:g
    % set movie boolean
    if i == movieindex
        movie = 1;
    else
        movie = 0;
    end
    
    xs = 0; % minimum x coordinate
    Lx = (1.5+0.75+1.75); % maximum x coordinate
    ys = 0; % minimum y coordinate
    Ly = (1.0+0.75+0.25); % maximum y coordinate
    
    h = Lx/M(i); % define step size
    
    % define x and y coordinates for u staggered mesh
    xf = linspace(xs, Lx, M(i)+1);
    yc = linspace(ys-h/2, Ly+h/2, N(i)+2);
    % define x and y coordinates for v staggered mesh
    xc = linspace(xs-h/2, Lx+h/2, M(i)+2);
    yf = linspace(ys, Ly, N(i)+1);
    % define x and y coordinates for cell centered mesh with 3 ghost cells
    xc3 = linspace(xs-5/2*h, Lx+5/2*h, M(i)+6);
    yc3 = linspace(ys-5/2*h, Ly+5/2*h, N(i)+6);
    
    Re = 20; % Reynolds number
    Sc = 2; % Schmidt number
    
    if movieindex == 0
        a = 5;
        b = 15;
        tout = [a,b];
    else
        tout = [1e-6 1, 4, 5, 5.5, 7.5, 10, 12, 15]; % time units to plot solution
    end
    Lt = length(tout); % define initial length of tout
    t = 0; % initialize time to zero
    tf = 15; % set final time
    
    % the fluid in the chamber is initially at rest
    u = zeros(M(i)+1, N(i)+2);
    v = zeros(M(i)+2, N(i)+1);
    % with Y(x, y, t=0) = 0
    Y = zeros(M(i)+6, N(i)+6);
    
    % preallocate hyperbolic terms
    Hu = zeros(size(u));
    Hv = zeros(size(v));
    Hy = zeros(size(Y));
    
    % preallocate Lagrangian operator
    phi = zeros(M(i)+2, N(i)+2);
    
    % source terms set to zero
    Qu = zeros(size(u)); Qv = zeros(size(v)); QY = zeros(size(Y));
    
    % preallocate time,deltat,k,S,
    time = zeros(10000, 1);
    deltat = time;
    k = time;
    S = time;
    
    if movie == 1
        % open movie files
        open(umovie);
        open(vmovie);
        open(Ymovie);
        open(Rmovie);
    end
    
    % define CFL
    CFL = 0.85;
    % counter for while loop
    count = 1;
    % counter for number of iterations
    iter = 1;
    % maximum number of V-cycle iterations
    n = 5e3;
    % error tolerance
    epsilon = 1e-10;
    
    % loop until final time is reached
    while t < tf
        % display current time to monitor runtime
        fprintf("%x\t%x\n",i,t);
        % add time to time vector
        time(iter,1) = t;
        % store t^(n-1) values
        u0 = u; v0 = v; Hu0 = Hu; Hv0 = Hv; Hy0 = Hy;
        % calculate stable time step
        [dt, outputFlag] = calcDt561(t,tout(count),u0,v0);
        deltat(iter,1) = dt;
        % calculate solution for u, v, Y
        [Hu, Hv] = hyperbolic_uv_2D(u0,v0);
        Hu = 3/2*Hu-1/2*Hu0;
        Hv = 3/2*Hv-1/2*Hv0;
        Hy = hyperbolic_Y_WENO_2D(Y,u0,v0,u,v,dt);
        % CN1 and apply boundaries coniditions
        u = parabolic_CN1_2D_u(u,Qu+Hu,dt);
        v = parabolic_CN1_2D_v(v,Qv+Hv,dt);
        Y = parabolic_CN1_2D_Y3(Y,QY+Hy,dt);
        % apply boundary conditions to u, v, Y, and R
        u = bc_u(u, t);
        v = bc_v(v, t);
        Y = bc_Y3(Y, t);
        % CN2 and apply boundary coinditions
        u = parabolic_CN2_2D_u(u,Qu+Hu,dt);
        v = parabolic_CN2_2D_v(v,Qv+Hv,dt);
        Y = parabolic_CN2_2D_Y3(Y,QY+Hy,dt);
        % apply boundary conditions to u, v, Y, and R
        u = bc_u(u,t);
        v = bc_v(v,t);
        Y = bc_Y3(Y,t);
        % correct outlet velocities
        [u,v] = correctOutlet(u,v);
        % define right hand side of Poisson equation
        f = 1/dt*calcDivV(u,v);
        % calculate right hand side of the Poisson equation
        phi = myPoisson(phi,f,h,n,epsilon);
        % project velocities
        [u, v] = projectV(u,v,phi,dt);
        % add value to S and k vectors
        S(iter, 1) = calcS3(Y);
        k(iter, 1) = calcK(u,v);
        % increase time
        t = t+dt;
        % increase iteration counter
        iter = iter+1;
        % plot at output time
        if outputFlag == 1
            if i == movieindex
                % plot u
                figure;
                pcolor(yc(2:N(i)+1),xf(2:M(i)),u(2:M(i), 2:N(i)+1));
                colormap(jet);
                shading interp;
                ylim([xs, Lx]);
                xlim([ys, Ly]);
                caxis([-1.5, 1.5]);
                colorbar;
                title("u Solution At t = "+num2str(tout(count))+" Seconds");
                xlabel("y");
                ylabel("x");
                view([90, 90]);
                set(gca, "XDir", "reverse")
                exportgraphics(gcf,"u"+num2str(tout(count)+".png"),"Resolution",1000);
                % plot v
                figure;
                pcolor(yf(2:N(i)),xc(2:M(i)+1),v(2:M(i)+1, 2:N(i)));
                colormap(jet);
                shading interp;
                ylim([xs, Lx]);
                xlim([ys, Ly]);
                caxis([-1.5, 1.5]);
                colorbar;
                title("v Solution At t = "+num2str(tout(count))+" Seconds");
                xlabel("y");
                ylabel("x");
                view([90, 90]);
                set(gca, "XDir", "reverse")
                exportgraphics(gcf,"v"+num2str(tout(count)+".png"),"Resolution",1000);
                % plot Y
                figure;
                pcolor(yc3(4:N(i)+3),xc3(4:M(i)+3),Y(4:M(i)+3, 4:N(i)+3));
                colormap(hot);
                shading interp;
                ylim([xs, Lx]);
                xlim([ys, Ly]);
                caxis([0, 1]);
                colorbar;
                title("Y Solution At t = "+num2str(tout(count))+" Seconds");
                xlabel("y");
                ylabel("x");
                view([90, 90]);
                set(gca, "XDir", "reverse")
                exportgraphics(gcf,"Y"+num2str(tout(count)+".png"),"Resolution",1000);
                % plot R
                figure;
                R = Y.*(1-Y);
                pcolor(yc3(4:N(i)+3),xc3(4:M(i)+3),R(4:M(i)+3, 4:N(i)+3));
                colormap(hot);
                shading interp;
                caxis([0, 0.25]);
                colorbar;
                title("R Solution At t = "+num2str(tout(count))+" Seconds");
                xlabel("y");
                ylabel("x");
                view([90, 90]);
                set(gca, "XDir", "reverse")
                exportgraphics(gcf,"R"+num2str(tout(count)+".png"),"Resolution",1000);
                % close all plots
                close all;
            elseif movieindex == 0
                T(i) = T(i)+deltat(iter-2,1)/2*S(iter-1, 1);
                KE(i) = KE(i)+deltat(iter-2,1)/2*k(iter-1, 1);
            end
            % increase counter
            count = count+1;
        elseif t > 5 && movieindex == 0
            T(i) = T(i)+dt*S(iter-1, 1);
            KE(i) = KE(i)+dt*k(iter-1, 1);
        end
        if movie == 1 && mod(iter-1,30) == 0
            % get frames for movies
            % plot u
            figure;
            pcolor(yc(2:N(i)+1),xf(2:M(i)),u(2:M(i), 2:N(i)+1));
            colormap(jet);
            shading interp;
            ylim([xs, Lx]);
            xlim([ys, Ly]);
            caxis([-1.5, 1.5]);
            colorbar;
            title("u Solution At t = "+num2str(time(iter-1,1))+" Seconds");
            xlabel("y");
            ylabel("x");
            view([90, 90]);
            set(gca, "XDir", "reverse")
            writeVideo(umovie,getframe(gcf));
            % plot v
            figure;
            pcolor(yf(2:N(i)),xc(2:M(i)+1),v(2:M(i)+1, 2:N(i)));
            colormap(jet);
            shading interp;
            ylim([xs, Lx]);
            xlim([ys, Ly]);
            caxis([-1.5, 1.5]);
            colorbar;
            title("v Solution At t = "+num2str(time(iter-1,1))+" Seconds");
            xlabel("y");
            ylabel("x");
            view([90, 90]);
            set(gca, "XDir", "reverse")
            writeVideo(vmovie,getframe(gcf));
            % plot Y
            figure;
            pcolor(yc3(4:N(i)+3),xc3(4:M(i)+3),Y(4:M(i)+3, 4:N(i)+3));
            colormap(hot);
            shading interp;
            ylim([xs, Lx]);
            xlim([ys, Ly]);
            caxis([0, 1]);
            colorbar;
            title("Y Solution At t = "+num2str(time(iter-1,1))+" Seconds");
            xlabel("y");
            ylabel("x");
            view([90, 90]);
            set(gca, "XDir", "reverse")
            writeVideo(Ymovie,getframe(gcf));
            % plot R
            figure;
            R = Y.*(1-Y);
            pcolor(yc3(4:N(i)+3),xc3(4:M(i)+3),R(4:M(i)+3, 4:N(i)+3));
            colormap(hot);
            shading interp;
            caxis([0, 0.25]);
            colorbar;
            title("R Solution At t = "+num2str(time(iter-1,1))+" Seconds");
            xlabel("y");
            ylabel("x");
            view([90, 90]);
            set(gca, "XDir", "reverse")
            writeVideo(Rmovie, getframe(gcf));
            % close all plots
            close all;
        end
    end
    
    
    if movie == 1
        % close movie files
        close(umovie);
        close(vmovie);
        close(Ymovie);
        close(Rmovie);
    end
    
    % plot S and k as functions of time
    if i == movieindex
        figure
        plot(time(1:iter-1), S(1:iter-1));
        title("S Versus Time");
        xlim([time(1) time(iter-1)]);
        xlabel("t");
        ylabel("S");
        exportgraphics(gcf,"Svt.png","Resolution",1000);
        figure
        plot(time(1:iter-1), k(1:iter-1));
        title("k Versus Time");
        xlim([time(1) time(iter-1)]);
        xlabel("t");
        ylabel("k");
        exportgraphics(gcf,"kvt.png","Resolution",1000);
        close all;
    end
end

% calculate T,KE
T = 1/10.*T;
KE = 1/10.*KE;