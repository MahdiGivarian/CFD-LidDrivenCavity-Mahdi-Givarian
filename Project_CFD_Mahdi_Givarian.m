% In The Name of God
% Course Name: Computational Fluid Dynamics (CFD)
% Instructor: Dr.Hamid Naderan
% By: Mahdi Givarian
% Student Number: 401129139

clc
clear
close all;

%%% Input of Viscosity
nu = input('Enter the Viscosity:  '); % Kinematic Viscosity

%%% GIVENS
N = [2064 1032 516 258 129]; L = 1; Wall_Velocity = 1; % Nodes; Domain Size; Velocity
rho = 1; mu = rho*nu; Re = (Wall_Velocity*L)/nu; % Density; Dynamic Viscosity; Reynolds number
dt = 0.000001; maxIt = 200000; maxe = 1e-7; % Time Step; Max iter; Max error
Delta = zeros(length(N)-1); successive_error = zeros(length(N)-1);
counter = 1;

%%% SETUP 2D GRID
for n = 1:length(N)
    Nx = N(n); Ny = Nx; h = L/(Nx-1);
    x = 0:h:L; y = 0:h:L;
    im = 1:Nx-2; i = 2:Nx-1; ip = 3:Nx; jm = 1:Ny-2; j = 2:Ny-1; jp = 3:Ny;
    
    %%% PRELOCATE MATRIXES
    Vo = zeros(Nx,Ny); St = Vo; Vop = Vo; u = Vo; v = Vo;
    
    %%% SOLVE LOOP SIMILAR TO GAUSS-SIEDEL METHOD
    for iter1 = 1:maxIt
        %%% CREATE BOUNDARY CONDITIONS
        Vo(1:Nx,Ny) = -2*St(1:Nx,Ny-1)/(h^2) - 2*Wall_Velocity/h; % Top
        Vo(1:Nx,1)  = -2*St(1:Nx,2)/(h^2);                        % Bottom
        Vo(1,1:Ny)  = -2*St(2,1:Ny)/(h^2);                        % Left
        Vo(Nx,1:Ny) = -2*St(Nx-1,1:Ny)/(h^2);                     % Right
        
        %%% PARTIALLY SOLVE VORTICITY TRANSPORT EQUATION
        Vop = Vo;
        Vo(i,j) = Vop(i,j) + ...
        (-1*(St(i,jp)-St(i,jm))/(2*h) .* (Vop(ip,j)-Vop(im,j))/(2*h)+...
        (St(ip,j)-St(im,j))/(2*h) .* (Vop(i,jp)-Vop(i,jm))/(2*h)+...
        nu*(Vop(ip,j)+Vop(im,j)-4*Vop(i,j)+Vop(i,jp)+Vop(i,jm))/(h^2))*dt;
        
        %%% PARTIALLY SOLVE ELLIPTICAL VORTICITY EQUATION FOR STREAM FUNCTION
        St(i,j) = (Vo(i,j)*h^2 + St(ip,j) + St(im,j) + St(i,jp) + St(i,jm))/4;
        
        %%% CHECK FOR CONVERGENCE (Steady State)
        if iter1 > 10
        error1 = max(max(Vo - Vop));
        if iter1 > 1000
        Error_residual(counter,n) = error1;
        counter = counter+1;
        end
        if error1 < maxe
        break;
        end
        end
    end
    
    %%% CREATE VELOCITY FROM STREAM FUNCTION
    u(2:Nx-1,Ny) = Wall_Velocity;
    u(i,j) = (St(i,jp)-St(i,jm))/(2*h); v(i,j) = (-St(ip,j)+St(im,j))/(2*h);
    
    % SOLVE POISSON EQUATION FOR PRESSURE
    P = ones(Nx,Ny);     % Initial Guess
    for iter2 = 1:maxIt
        %%% CREATE BOUNDARY CONDITIONS for Interior and Corner Points
        P(2:Nx,Ny)   = P(1:Nx-1,Ny)-mu/2*(3*Vo(2:Nx,Ny)-4*Vo(2:Nx,Ny-1)+Vo(2:Nx,Ny-2));        % Top
        P(2:Nx,1)    = P(1:Nx-1,1)-mu/2*(-3*Vo(2:Nx,1)+4*Vo(2:Nx,2)-Vo(2:Nx,3));               % Bottom
        P(1,2:Ny)    = P(1,1:Ny-1)+mu/2*(-3*Vo(1,2:Ny)+4*Vo(2,2:Ny)-Vo(3,2:Ny));               % Left
        P(Nx,2:Ny-1) = P(Nx,1:Ny-2)+mu/2*(3*Vo(Nx,2:Ny-1)-4*Vo(Nx-1,2:Ny-1)+Vo(Nx-2,2:Ny-1));  % Right
        
        %%% PARTIALLY SOLVE POISSON PRESSURE EQUATION
        Pp = P;
        for i = 2:Nx-1
            for j = 2:Ny-1
                P(i,j) = (rho*(-2*(St(i+1,j)+St(i-1,j)-2*St(i,j))*(St(i,j+1)+St(i,j-1)-2*St(i,j))/(h^2)+((St(i+1,j+1)-St(i-1,j+1)-St(i+1,j-1)+St(i-1,j-1))^2)/(8*h^2))+P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1))/4;
            end
        end
        
        %%% CHECK FOR CONVERGENCE
        if iter2 > 10
        error2 = max(max(P - Pp));
        if error2 < maxe
        break;
        end
        end
    end
    
    if n==1
        u_1 = u;
    elseif n==2
        u_2 = u;
        n_2 = Nx;
        Delta(n-1) = h;
        % Calculating the Error
        sum = 0;
        for i = 1:n_2
            for j = 1:n_2
                Error = abs(u_2(i,j)-u_1(2*i,2*j));
                % Calculating the Norm2 of the error L_2(E)
                sum = sum+Error^2;
            end
        end
        successive_error(n-1) = sqrt(sum/n_2);
    
    elseif n==3
        u_3 = u;
        n_3 = Nx;
        Delta(n-1) = h;
        % Calculating the Error
        sum = 0;
        for i = 1:n_3
            for j = 1:n_3
                Error = abs(u_3(i,j)-u_2(2*i,2*j));
                % Calculating the Norm2 of the error L_2(E)
                sum = sum+Error^2;
            end
        end
        successive_error(n-1) = sqrt(sum/n_3);
        
    elseif n==4
        u_4 = u;
        n_4 = Nx;
        Delta(n-1) = h;
        % Calculating the Error
        sum = 0;
        for i = 1:n_4
            for j = 1:n_4
                Error = abs(u_4(i,j)-u_3(2*i,2*j));
                % Calculating the Norm2 of the error L_2(E)
                sum = sum+Error^2;
            end
        end
        successive_error(n-1) = sqrt(sum/n_4);
        
    elseif n==5
        u_5 = u;
        n_5 = Nx;
        Delta(n-1) = h;
        % Calculating the Error
        sum = 0;
        for i = 1:n_5
            for j = 1:n_5
                Error = abs(u_5(i,j)-u_4(2*i,2*j));
                % Calculating the Norm2 of the error L_2(E)
                sum = sum+Error^2;
            end
        end
        successive_error(n-1) = sqrt(sum/n_5);    
    end

end

Iteration = 1001:iter1;
[X,Y] = meshgrid(x,y);
M = 1000; xstart = max(x)*rand(M,1); ystart = max(y)*rand(M,1);

%%% PLOTS
cm = hsv(ceil(100/0.7)); 
cm = flipud(cm(1:100,:));
figure(1) 
contourf(x,y,u',23,'LineColor','none')
title('U-velocity')
xlabel('x-location')
ylabel('y-location')
axis('equal',[0 L 0 L]) 
colormap(cm) 
colorbar('westoutside')

figure(2) 
contourf(x,y,v',23,'LineColor','none')
title('V-velocity')
xlabel('x-location') 
ylabel('y-location')
axis('equal',[0 L 0 L]) 
colormap(cm) 
colorbar('westoutside')

figure(3) 
plot(y,u(round(Nx/2),:))
title('Centerline x-direction velocity')
xlabel('y/L') 
ylabel('u/U') 
axis('square') 
xlim([0 L]) 
grid on

figure(4) 
h = streamline(X,Y,u',v',xstart,ystart,[0.1, 200]);
title('Stream Function') 
xlabel('x-location') 
ylabel('y-location')
axis('equal',[0 L 0 L]) 
set(h,'color','k')
hold on

figure(5)
contourf(X,Y,St',40)
view(2)
shading interp
xlabel('x')
ylabel('y')
title(['Stream Function for Re=',num2str(Re)])
colorbar
grid on

figure(6)
contourf(X,Y,Vo',240)
view(2)
shading interp
xlabel('x')
ylabel('y')
title(['Vorticity for Re=',num2str(Re)])
colorbar
grid on

figure(7)
contourf(x,y,P',30,'LineColor','none')
view(2)
shading interp
xlabel('x')
ylabel('y')
title(['Pressure for Re=',num2str(Re)])
colorbar
grid on

figure(8)
loglog(Delta,successive_error)
xlabel('h')
ylabel('Successive Error')

figure(9)
plot(Iteration,Error_residual(:,5))
xlabel('Iterations')
ylabel('Residual Error')
title(['Residual Error (1) for Re=',num2str(Re)])
grid on

figure(10)
plot(Iteration(24001:iter1),Error_residual(24001:iter1,5))
xlabel('Iterations')
ylabel('Residual Error')
title(['Residual Error (2) for Re=',num2str(Re)])
grid on
