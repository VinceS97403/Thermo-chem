% Thermal code 1D

clear all
close all
clf

% Define parameters
max_depth = 15e3; %km
k = 3; % W/m/K
Cp = 1000; % J/kg/K
rho = 3000; % kg.m^-3

T_back = 1000; %C
T_intr = 1300; %C

kappa = k/(rho*Cp);

% Intrusion

begin = 10e3; %km intrusion beginning
finish = 11e3; %km intrusion ending

% Define numerical grid space
ysize = max_depth;
numy = 100;
dy = ysize / (numy-1);

y = 0:dy:ysize;

% Set time stepping
dt = (dy^2)/(2*kappa);
%dt = 1e6;
numt = 1000; % number of nodes

% Initiate vector and matrix
R = zeros(numy,1);
L = zeros (numy,numy);

T0 = zeros(numy,1);
% Intrusion temperature profile
for i = 1:1:numy
    T0(i,1) = T_back;
    
    if (y(1,i)>begin && y(1,i)<finish)
        T0(i,1) = T_intr;
    else
        T0(i,1) = T_back;
    end
end
pro = T0;

%% Fill in R and L
% Boundary condition 1: T1 = 1000
R(1,1) = T0(1,1);
L(1,1) = 1;

% Boundary condition 2: T_end = 1000
R(numy,1) = T_back;
L(numy, numy) = 1;

% Loop over time until 3000 years;
timesum = 0;
T = T0;
figure(1)
while timesum < 3000*365*24*3600 %3000 years in seconds
       % Loop over space, as in poisson equation code
       for iy = 2:1:numy-1
           R(iy,1)= T0(iy,1)/dt;

           L(iy, iy-1) = -kappa/(dy^2);
           L(iy,iy) = 1/dt + 2*kappa/(dy^2);
           L(iy,iy+1) = -kappa/(dy^2);
       
       
       end 
       
       % Calculate solution
       T = L \ R;
       S = T;
       
       % Increment time
       timesum = dt + timesum;
       
      T0 = T;

      % Add solution to plot
%     Temp=plot(T, -y, 'LineWidth', 1);
%     drawnow
%    
%    if it == 1
%        Legend{it+1}=string(it)+' year';
%    else
%        Legend{it+1}=string(it)+' years';
%    end

%     if timesum > 3600*24*365*3000 %3000 years in seconds
%        timesumy = timesum /(3600*24*365);
%        Temp=plot(T, -y, 'LineWidth', 1);
%        drawnow
%        peak = max(T)
%     else
%     end
    plot(T, -y,'r', 'LineWidth', 1);
    xlim([900 1400])
    ylim([-15e3 0])
    xlabel('Temperature (C)')
    ylabel('Depth (km)')
    title(['Time= ' num2str(timesum/(3600*24*365)) +' years'])
    drawnow
 %end
end

hold on

Tempro=plot(pro,-y, 'b', 'LineWidth',1);
Legend{2}='Time of intrusion' ;

peak = max(T)

timesumy = timesum /(3600*24*365);

Legend{1}=string(timesumy)+' years';

legend(Legend);