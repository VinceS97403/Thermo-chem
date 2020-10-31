% 1D Subduction channel

clear all
clf

% Define parameters
L = 5e3; %m
eta = 1e23; %Pa.s
rho = 3000; %kg.m^-3
dPdY = 100000; %Pa.m^-1
gy = 9.81; %m.s^-2

% Set velocity boundary conditions
v_bc_up = -1e-12; %m.s^-1
v_bc_down = 0; %m.s^-1

% Calculate right hand side
RHS = (1/eta)*(dPdY - rho*gy);

% Add analytical solution for benchmarking
% x_ana = 0:0.001:L;
% vy_ana = (1/(2*eta))*(dPdY-rho*gy)*(x_ana.^2-L*x_ana);

% Define numerical grid in space
xsize = L; % Total length
numx = 20; % Number of nodes / Number of lines of the vector
dx = xsize / (numx-1); % Distance between the nodes

x = 0:dx:xsize; % Position of the nodes

% Initiate vector and matrix
% (row, column)
R = zeros(numx,1);
L = zeros(numx, numx);

%% Fill in R and L
% Boundary conditions 1: vy1 = 0
R(1,1) = v_bc_up;
L(1,1) = 1;

% Boundary condition N: vyN = 0
R(numx,1) = v_bc_down;
L(numx,numx) = 1; % Can also type L(end,end)

% Intermediate nodal points: coefficient in the stiffness matrix
for i=2:1:numx-1
   R(i,1)=RHS;
   
   L(i,i-1)=1/(dx^2);
   L(i,i)=-2/(dx^2);
   L(i,i+1)=1/(dx^2);
end

% Calculate solution
S = L \ R;

vy = -S;

figure(1);
Poisson = plot(x,vy, 'b.', 'MarkerSize', 40, 'LineWidth', 60);
xlabel('Relative position along the width of the channel (m)')
ylabel('Velocity (m/s)')
title('Mechanical model of rock flow pattern in a subduction channel', 'FontSize', 16)

hold on;

%Analytical = plot (x_ana,-vy_ana, 'r', 'LineWidth', 4);

hold on

Poisson_curve = plot(x,vy, 'k-','LineWidth', 2);

legend({'Poisson', 'Analytical'})