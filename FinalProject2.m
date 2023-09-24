clear all;
clc;
%N = [8 16 32];
%Re = [10 100 400];
N = 32;
Re = 400;
for l=1:length(N)
for r =1:length(Re)
%% Given
t_f = 3; % final time
dt = 0.01; % time step
u0 = 1; % velocity of lid
L_x = 1; % x-length of cavity
L_y = 1; % y-length of cavity
N_x = N(l); % # of x-nodes
N_y = N(l); % # of y-nodes
dx = L_x/N_x; % spatial-x step
dy = L_y/N_y; % spatial-y step
beta = dx/dy;
alpha = 1/(2*(1+beta^2));
w_PSOR = 0.9; % relaxation parameter for PSOR
ERROR_TOL = 0.05; % error tolerance
Re_r = Re(r); % Reynold's Number
%% Define Initial Conditions
u = zeros(N_x,N_y); % Define u velocity field
u(1, :) = 1; % initialize lid velocity
v = zeros(N_x,N_y); % Define v velocity field
% Initialize tridiagonal Matrix
tri_x = zeros(N_x,N_y);
tri_y = zeros(N_x,N_y);
% Initialize vorticity matrices
w_x = zeros(N_x, N_y); % X-sweep
w_y = zeros(N_x, N_y); % Y-sweep
w = zeros(N_x, N_y); % final vorticity matrix
% Initialize stream-function matrix
psi = zeros(N_x, N_y);
%% Boundary Conditions
w_x(:,1) = 2*(psi(:,1)-psi(:,2))/dy^2;
w_x(:,end) = 2*(psi(:,end)-psi(:,end-1))/dx^2;
w_x(end,:) = 2*(psi(end,:)-psi(end-1,:))/dx^2;
w_x(1,:) = 2*(psi(1,:)-psi(2,:))/dy^2 - 2*u0/dy;

w_y(:,1) = 2*(psi(:,1)-psi(:,2))/dy^2;
w_y(:,end) = 2*(psi(:,end)-psi(:,end-1))/dx^2;
w_y(end,:) = 2*(psi(end,:)-psi(end-1,:))/dx^2;
w_y(1,:) = 2*(psi(1,:)-psi(2,:))/dy^2 - 2*u0/dy;

w(:,1) = 2*(psi(:,1)-psi(:,2))/dy^2;
w(:,end) = 2*(psi(:,end)-psi(:,end-1))/dx^2;
w(end,:) = 2*(psi(end,:)-psi(end-1,:))/dx^2;
w(1,:) = 2*(psi(1,:)-psi(2,:))/dy^2 - 2*u0/dy;

% Define coefficents for tridiagonal matrix
d_x = 1/Re_r*dt/dx^2;
d_y = 1/Re_r*dt/dy^2;
B_x = 1 + d_x;
B_y = 1 + d_y;

n = 1; % iteration level
t = 0; % current time
current_error = 1;
while t < t_f%current_error > ERROR_TOL % loop until error is within the tolerance
    %% Boundary Conditions at next iteration level
    % Stream-Function
    psi(:,1, n+1) = psi(:,1, n); 
    psi(1,:, n+1) = psi(1,:, n); 
    psi(:,end, n+1) = psi(:,end, n); 
    psi(end,:, n+1) = psi(end,:, n);
    
    % Velocity Fields
    u(1,:, n+1) = u(1,:,n);
    u(end,:, n+1) = u(end,:,n);
    u(:,1, n+1) = u(:,1,n);
    u(:,end, n+1) = u(:,end,n);
    
    v(1,:, n+1) = v(1,:,n);
    v(end,:, n+1) = v(end,:,n);
    v(:,1, n+1) = v(:,1,n);
    v(:,end, n+1) = v(:,end,n);
    
    % Vorticity Matrices
    w_x(:,1,n+1) = 2*(psi(:,1,n+1)-psi(:,2,n+1))/dy^2;
    w_x(:,end,n+1) = 2*(psi(:,end,n+1)-psi(:,end-1,n+1))/dx^2;
    w_x(end,:,n+1) = 2*(psi(end,:,n+1)-psi(end-1,:,n+1))/dx^2;
    w_x(1,:,n+1) = 2*(psi(1,:,n+1)-psi(2,:,n+1))/dy^2 - 2*u0/dy;
    
    w_y(:,1,n+1) = 2*(psi(:,1,n+1)-psi(:,2,n+1))/dy^2;
    w_y(:,end,n+1) = 2*(psi(:,end,n+1)-psi(:,end-1,n+1))/dx^2;
    w_y(end,:,n+1) = 2*(psi(end,:,n+1)-psi(end-1,:,n+1))/dx^2;
    w_y(1,:,n+1) = 2*(psi(1,:,n+1)-psi(2,:,n+1))/dy^2 - 2*u0/dy;

    w(:,1,n+1) = 2*(psi(:,1,n+1)-psi(:,2,n+1))/dy^2;
    w(:,end,n+1) = 2*(psi(:,end,n+1)-psi(:,end-1,n+1))/dx^2;
    w(end,:,n+1) = 2*(psi(end,:,n+1)-psi(end-1,:,n+1))/dx^2;
    w(1,:,n+1) = 2*(psi(1,:,n+1)-psi(2,:,n+1))/dy^2 - 2*u0/dy;

    %% Solve Vorticity-Transport Equation by ADI Method
    % 1st define tridiagonal matrices' elements
    for i = 1:N_y
        for j = 1:N_x
            c_x(i,j) = u(i,j,n)*dt/dx;
            A_x(i,j) = -1/2*(1/2*c_x(i,j) + d_x);
            C_x(i,j) = 1/2*(1/2*c_x(i,j) - d_x);
            c_y(i,j) = v(i,j,n)*dt/dy;
            A_y(i,j) = -1/2*(1/2*c_y(i,j) + d_y);
            C_y(i,j) = 1/2*(1/2*c_y(i,j) - d_y);
            if i == j + 1 % Lower Diagonal
                tri_x(i,j) = A_x(i,j);
                tri_y(i,j) = A_y(i,j);
            elseif i == j % Main Diagonal
                tri_x(i,j) = B_x;
                tri_y(i,j) = B_y;
            elseif i == j - 1 % Upper Diagonal
                tri_x(i,j) = C_x(i,j);
                tri_y(i,j) = C_y(i,j);
            end
        end
    end
    % Perform X-SWEEP
    for i = 2:N_y-1 % Iterate bottom to top
        D_x = [];
        for j = 2:N_x-1 % Iterate left to right
            % initialize D matrix of "Ax = D"
            if j == 2
                D_ij = (c_y(i,j)/2+d_y)/2 * w(i,j-1,n) + (1-d_y)*w(i,j,n) + (-c_y(i,j)/2+d_y)/2*w(i,j+1,n) - A_x(i,j)*w(i,j-1,n);
            elseif j == N_x-1
                D_ij = (c_y(i,j)/2+d_y)/2 * w(i,j-1,n) + (1-d_y)*w(i,j,n) + (-c_y(i,j)/2+d_y)/2*w(i,j+1,n) - C_x(i,j)*w(i,j+1,n);
            else
                D_ij = (c_y(i,j)/2+d_y)/2 * w(i,j-1,n) + (1-d_y)*w(i,j,n) + (-c_y(i,j)/2+d_y)/2*w(i,j+1,n);
            end
            D_x = [D_x D_ij];
        end
        % Solve for vorticity at ith row using Matlab '/' matrix operator
        w_i_x = D_x/tri_x(2:end-1, 2:end-1);
        w_x(i, 2:end-1, n+1) = w_i_x;
    end
    % Y-SWEEP
    % 1st define tridiagonal matrices' elements
    % Perform Y-SWEEP
    for j= 2:N_x-1 % Iterate left to right
        D_y = [];
        for i = 2:N_y-1 % Iterate bottom to top
            % initialize D matrix of "Ax = D"
            if i == 2
                D_ij = (c_x(i,j)/2+d_x)/2 * w_x(i-1,j,n) + (1-d_x)*w_x(i,j,n) + (-c_x(i,j)/2+d_x)/2*w_x(i+1,j,n) - A_y(i,j)*w_x(i-1,j,n+1);
            elseif i == N_y-1
                D_ij = (c_x(i,j)/2+d_x)/2 * w_x(i-1,j,n) + (1-d_x)*w_x(i,j,n) + (-c_x(i,j)/2+d_x)/2*w_x(i+1,j,n) - C_y(i,j)*w_x(i+1,j,n+1);
            else
                D_ij = (c_x(i,j)/2+d_x)/2 * w_x(i-1,j,n) + (1-d_x)*w_x(i,j,n) + (-c_x(i,j)/2+d_x)/2*w_x(i+1,j,n);
            end
            D_y = [D_y D_ij];
        end
        % Solve for vorticity at jth column using Matlab '/' matrix operator
        w_j_y = D_y/tri_y(2:end-1, 2:end-1);
        w_y(2:end-1, j, n+1) = w_j_y;
        w = w_y;
    end
    %% Solve Stream Function by PSOR
    PSOR_error = 1; % initialize current PSOR error
    count = n;
    while PSOR_error > ERROR_TOL % iterate until PSOR error is within the tolerance
        psi_old = psi(:,:,end-1);
        for i = 2:N_y-1
            for j = 2:N_x-1
                psi(i,j,n+1) = (1-w_PSOR)*psi(i,j,n) + w_PSOR*alpha*(w(i,j,n+1)*dx^2+ psi(i+1,j,n) + psi(i-1,j,n+1)+beta^2*(psi(i,j+1,n)+psi(i,j-1,n+1)));
            end
        end
        psi_old = psi(:,:,end-1);
        PSOR_error = max(max(abs(psi(:,:,end) - psi_old)));
    end
    n = count;
    %% Initialize new velocity fields
    for i = 2:N_y-1
        for j = 2:N_x-1
            u(i,j,n+1) = (psi(i,j+1,n)-psi(i,j-1,n))/(2*dy);
            v(i,j,n+1) = -(psi(i+1,j,n)-psi(i-1,j,n))/(2*dx);
        end
    end
    %% Update Counters and Error to Check for Convergence
    current_error = max(max(abs(w_y(2:end-1,2:end-1,n+1) - w_y(2:end-1,2:end-1,n))));
    n = n + 1; % increment iteration level
    t = t + dt; % increment time
end
w_flipped = flipud(w(:,:,end));
psi_flipped = flipud(psi(:,:,end));
u_flipped = flipud(u(:,:,end));
v_flipped = flipud(v(:,:,end));


%% Plots
Re_str = num2str(Re(r));
N_str = num2str(N(l));
figure(3*(l-1)+r) % Vorticity Plot
contourf(w_flipped)
hold on
%title(['Vorticity for Re = ', Re_str, ' and N = ', N_str, 'x', N_str]);
%figure(3*(l-1)+r+1) % Stream-function Plot
contour(psi_flipped)
hold on 
%title(['Stream Function for Re = ', Re_str, ' and N = ', N_str, 'x', N_str]);
%figure(3*(l-1)+r+2) % Velocity-field Plot
quiver(u_flipped(2:end-1,2:end-1), v_flipped(2:end-1,2:end-1), 1)
title(['Vorticity, Stream, and Velocity Fields for Re = ', Re_str, ' and N = ', N_str, 'x', N_str]);
end
end